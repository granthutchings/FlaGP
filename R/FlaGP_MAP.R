#' @title FlaGP MAP calibration
#'
#' @description MAP estimate of calibration parameters
#' @param flagp x an \code{flagp} object.
#' @param n.restarts size of random sample used to initialize the optimizer
#' @param init matrix of locations to initialize the optimizer
#' @param seed seed for reproducable initial conditions for the optimizer
#' @param opts list of options for crs::snomadr
#' @details Returns FlaGP data object
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
map = function(flagp,n.restarts=1,init=NULL,seed=1,
               opts=list("MAX_BB_EVAL" = 10000, "INITIAL_MESH_SIZE" = .1,"MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0),
               end.eta=50,delta.method='newGP',start.delta=6,end.delta=50,
               theta.prior='beta',theta.prior.params=c(2,2),
               ssq.prior='hcauchy',ssq.prior.params=c(.1),
               method='snomadr',make.cluster=F){

  p.t = flagp$XT.data$p.t
  if(!is.null(init)){
    # if init passed, check for correct dimensions
    if(ncol(init) != p.t+1){
      stop("init should be a matrix with p.t+1 columns")
    }
    n.restarts = nrow(init)
  } else{
    # if init not passed generate with LHS
    if(n.restarts>1){
      set.seed(seed)
      init = rbind(init,lhs::maximinLHS(max(2,n.restarts),p.t+1))
      init[,p.t+1] = init[,p.t+1]*.1
    } else{
      init = matrix(c(rep(.5,p.t),.01),nrow=1)
    }
  }
  time = proc.time()

  if(method=='snomadr'){
    out = NULL
    if(n.restarts>1){
      if(make.cluster){
        cores = min(n.restarts,parallel::detectCores()-1)
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl,cores)
      }
      out = foreach::foreach(i=1:n.restarts) %dopar% crs::snomadr(fit_model_map, n = p.t+1,
                                                                  bbin = rep(0, p.t+1), bbout = 0,
                                                                  x0 = init[i,],
                                                                  lb = rep(0, p.t+1),
                                                                  ub = rep(1, p.t+1),
                                                                  random.seed = seed,
                                                                  opts = opts,
                                                                  print.output=F,
                                                                  flagp = flagp,
                                                                  end.eta=end.eta,
                                                                  delta.method=delta.method,
                                                                  start.delta=start.delta,
                                                                  end.delta=end.delta,
                                                                  theta.prior=theta.prior,
                                                                  theta.prior.params=theta.prior.params,
                                                                  ssq.prior=ssq.prior,
                                                                  ssq.prior.params=ssq.prior.params)
    } else{
      out = crs::snomadr(fit_model_map, n = p.t+1, bbin = rep(0, p.t+1), bbout = 0,
                         x0 = init,
                         lb = rep(0, p.t+1),
                         ub = rep(1, p.t+1),
                         random.seed = seed,
                         opts = opts,
                         print.output=F,
                         flagp = flagp,
                         end.eta=end.eta,
                         delta.method=delta.method,
                         start.delta=start.delta,
                         end.delta=end.delta,
                         theta.prior=theta.prior,
                         theta.prior.params=theta.prior.params,
                         ssq.prior=ssq.prior,
                         ssq.prior.params=ssq.prior.params)
    }

    if(n.restarts>1){
      solutions = matrix(nrow=n.restarts,ncol=p.t+1)
      objectives = rep(0,n.restarts)
      for(i in 1:n.restarts){
        solutions[i,] = out[[i]]$solution
        objectives[i] = out[[i]]$objective
      }
      theta.hat = solutions[which.min(objectives),1:p.t]
      ssq.hat = solutions[which.min(objectives),p.t+1]
    } else{
      solutions = out$solution
      theta.hat = out$solution[1:p.t]
      ssq.hat = out$solution[p.t+1]
    }
  } else{
    # regular optim
    if(n.restarts>1){
      if(make.cluster){
        cores = min(n.restarts,parallel::detectCores()-1)
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl,cores)
      }
      out = foreach::foreach(i=1:n.restarts) %dopar% optim(par = init[i,],
                                                           fn = fit_model_map,
                                                           # lower = rep(0, p.t+1),
                                                           # upper = c(rep(1-1e8, p.t),Inf),
                                                           # method = 'L-BFGS-B',
                                                           flagp = flagp,
                                                           end.eta=end.eta,
                                                           delta.method=delta.method,
                                                           start.delta=start.delta,
                                                           end.delta=end.delta,
                                                           theta.prior=theta.prior,
                                                           theta.prior.params=theta.prior.params,
                                                           ssq.prior=ssq.prior,
                                                           ssq.prior.params=ssq.prior.params,
                                                           hessian = T)
    } else{
      out = optim(par = init,
                   fn = fit_model_map,
                   #lower = rep(0, p.t+1),
                   #upper = c(rep(1-1e8, p.t),Inf),
                   #method = 'L-BFGS-B',
                   flagp = flagp,
                   end.eta=end.eta,
                   delta.method=delta.method,
                   start.delta=start.delta,
                   end.delta=end.delta,
                   theta.prior=theta.prior,
                   theta.prior.params=theta.prior.params,
                   ssq.prior=ssq.prior,
                   ssq.prior.params=ssq.prior.params,
                   hessian = T)
    }
    if(n.restarts>1){
      solutions = matrix(nrow=n.restarts,ncol=p.t+1)
      objectives = rep(0,n.restarts)
      for(i in 1:n.restarts){
        solutions[i,] = out[[i]]$par
        objectives[i] = out[[i]]$value
      }
      theta.hat = solutions[which.min(objectives),1:p.t]
      ssq.hat = solutions[which.min(objectives),p.t+1]
      Cov = chol2inv(chol(out[[which.min(objectives)]]$hessian))
    } else{
      solutions = out$par
      theta.hat = out$par[1:p.t]
      ssq.hat = out$par[p.t+1]
      Cov = chol2inv(chol(out$hessian)) # not -H because we did minimzation of negll not maximination of ll
    }
  }
  opt.time = proc.time() - time

  if(make.cluster){
    parallel::stopCluster(cl)
  }
  return = fit_model(theta.hat,ssq.hat,flagp,lite=F,sample=F,end.eta=end.eta,
                     delta.method=delta.method,start.delta=start.delta,end.delta=end.delta,
                     theta.prior=theta.prior,theta.prior.params=theta.prior.params,
                     ssq.prior=ssq.prior,ssq.prior.params=ssq.prior.params,map=T)
  if(method=='optim')
    return$Cov = Cov
  return$delta$method = delta.method
  return$opt = out
  return$theta.hat = theta.hat; return$theta = NULL
  return$ssq.hat = ssq.hat; return$ssq = NULL
  return$solutions = solutions
  return$time = opt.time
  return$init = init[1:n.restarts,]
  class(return) = c('map',class(return))
  return(return)
}
