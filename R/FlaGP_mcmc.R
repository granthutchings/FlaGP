mcmc_mh = function(flagp,
                   t.init=rep(.5,flagp$XT.data$p.t),
                   ssq.init=.01,
                   prop.cov=diag((.5/3)^2,flagp$XT.data$p.t+1),
                   n.samples=10000,n.burn=1000,thin=1,
                   adapt.par = c(100,50,.5,1000),
                   end.eta=50,
                   delta.method='newGP',start.delta=6,end.delta=50,
                   theta.prior='beta',theta.prior.params=c(2,2),
                   ssq.prior='hcauchy',ssq.prior.params=c(.1),
                   prev.samples=NULL,
                   recalc.llh=F,
                   sample=as.logical(ifelse(flagp$bias,F,T)),
                   verbose=T,loglogit=F,single_site=F){
  llh = function(par,full=F,opt=F,sample=T,...){
    t.curr = par[1:flagp$XT.data$p.t]
    ssq.curr = par[flagp$XT.data$p.t+1]
    if(ssq.curr<0 | any(t.curr<0) | any(t.curr>1)){
      if(full){
        return(list(ll=-Inf))
      } else{
        return(-Inf)
      }
    }
    if(opt){
      llt = fit_model(t.curr,ssq.curr,flagp,F,F,end.eta,delta.method,start.delta,end.delta,F,
                      theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
    } else{
      llt = fit_model(t.curr,ssq.curr,flagp,F,sample=sample,end.eta,delta.method,start.delta,end.delta,F,
                      theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
    }

    if(full){
      return(llt)
    } else{
      return(llt$ll)
    }
  }

  flagp$prior = list(delta.method=delta.method,start.delta=start.delta,end.delta=end.delta,
                     theta.prior=theta.prior,theta.prior.params=theta.prior.params,
                     ssq.prior=ssq.prior,ssq.prior.params=ssq.prior.params)
  if(loglogit){
    if(single_site){
      cat('Single-Site Metropolis on unconstrained (logit/log transformed) parameter space\n')
      mcmc = Metro_Hastings_Stochastic_Logit_OneAtATime(flagp,llh,c(t.init,ssq.init),prop.cov,NULL,n.samples,n.burn,thin,adapt.par,recalc.llh,sample,!verbose)
    } else{
      cat('Joint Metropolis on unconstrained (logit/log transformed) parameter space\n')
      mcmc = Metro_Hastings_Stochastic_Logit(flagp,llh,c(t.init,ssq.init),prop.cov,NULL,n.samples,n.burn,thin,adapt.par,recalc.llh,sample,!verbose)
    }
  } else{
    cat('Sampling in bounded space\n')
    mcmc = Metro_Hastings_Stochastic(flagp,llh,c(t.init,ssq.init),prop.cov,NULL,n.samples,n.burn,thin,adapt.par,recalc.llh,sample,!verbose)
  }

  # store some other stuff
  mcmc$dat = list(
    loglogit = loglogit,
    prior = flagp$prior,
    recalc.llh = recalc.llh,
    sample = sample,
    end.eta = end.eta
  )

  return(mcmc)
}

mcmc_mh_add = function(mcmc,flagp,verbose,n.samples=10000){
  llh = function(par,full=F,opt=F,sample=T,...){
    p.t = length(par) - 1
    t.curr = par[1:p.t]
    ssq.curr = par[p.t+1]
    if(ssq.curr<0 | any(t.curr<0) | any(t.curr>1)){
      if(full){
        return(list(ll=-Inf))
      } else{
        return(-Inf)
      }
    }
    if(opt){
      llt = fit_model(t.curr,ssq.curr,flagp,F,F,flagp$prior$end.eta,flagp$prior$delta.method,flagp$prior$start.delta,end.delta,F,
                      flagp$prior$theta.prior,flagp$prior$theta.prior.params,flagp$prior$ssq.prior,flagp$prior$ssq.prior.params)
    } else{
      llt = fit_model(t.curr,ssq.curr,flagp,F,sample,flagp$prior$end.eta,flagp$prior$delta.method,flagp$prior$start.delta,flagp$prior$end.delta,F,
                      flagp$prior$theta.prior,flagp$prior$theta.prior.params,flagp$prior$ssq.prior,flagp$prior$ssq.prior.params)
    }

    if(full){
      return(llt)
    } else{
      return(llt$ll)
    }
  }
  n.samp.prev = nrow(mcmc$t.samp)
  t.init = c(mcmc$t.samp[n.samp.prev,])
  ssq.init = mcmc$ssq.samp[n.samp.prev]
  prop.cov = mcmc$prop.cov
  recalc.llh = mcmc$dat$recalc.llh
  sample = mcmc$dat$sample
  loglogit = mcmc$dat$loglogit
  flagp$prior = mcmc$dat$prior
  flagp$prior$end.eta = mcmc$dat$end.eta
  n.burn = 0
  thin = 1
  adapt.par = c(0,0,0,0)

  if(loglogit){
    cat('Sampling unconstrained transformed space\n')
    mcmc_new = Metro_Hastings_Stochastic_Logit(flagp,llh,c(t.init,ssq.init),prop.cov,NULL,n.samples,n.burn,thin,adapt.par,recalc.llh,sample,!verbose)
  } else{
    cat('Sampling in bounded space\n')
    mcmc_new = Metro_Hastings_Stochastic(flagp,llh,c(t.init,ssq.init),prop.cov,NULL,n.samples,n.burn,thin,adapt.par,recalc.llh,sample,!verbose)
  }

  mcmc$t.samp = rbind(mcmc$t.samp,mcmc_new$t.samp)
  mcmc$ssq.samp = c(mcmc$ssq.samp,mcmc_new$ssq.samp)
  mcmc$ll.samp = c(mcmc$ll.samp,mcmc_new$ll.samp)
  mcmc$eta = c(mcmc$eta,mcmc_new$eta)
  mcmc$delta = c(mcmc$delta,mcmc_new$delta)
  mcmc$prop.cov = mcmc_new$prop.cov
  mcmc$time = mcmc$time + mcmc_new$time
  mcmc$acpt.ratio = ((mcmc$acpt.ratio * mcmc$n.samples) + (mcmc_new$acpt.ratio * mcmc_new$n.samples)) / (mcmc$n.samples + mcmc_new$n.samples)
  mcmc$n.samples = mcmc$n.samples + mcmc_new$n.samples
  mcmc$samp.all = rbind(mcmc$samp.all,mcmc_new$samp.all)
  return(mcmc)
}
# addapted MHadaptive::Metro_Hastings for stochastic llh and fast resampling of llh
Metro_Hastings_Stochastic = function (flagp, li_func, pars, prop_sigma = NULL, par_names = NULL,
                                      iterations = 50000, burn_in = 1000, thin = 1, adapt_par = c(100, 20, 0.5, 0.75), recalc = TRUE, sample = TRUE, quiet = FALSE, ...)
{
  ptm = proc.time()[3]
  if(is.null(sample))
    sample = as.logical(ifelse(flagp$bias,F,T))
  if (!is.finite(li_func(pars, sample=sample,...)))
    stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")
  if (is.null(par_names))
    par_names <- letters[1:length(pars)]
  if (!is.null(dim(prop_sigma))) {
    if ((dim(prop_sigma)[1] != length(pars) || dim(prop_sigma)[2] !=
         length(pars)) && !is.null(prop_sigma))
      stop("prop_sigma not of dimension length(pars) x length(pars)")
  }
  if (is.null(prop_sigma)) {
    if (length(pars) != 1) {
      fit <- optim(pars, li_func, control = list(fnscale = -1),
                   hessian = TRUE, opt=T, sample=sample,...)
      fisher_info <- solve(-fit$hessian)
      prop_sigma <- sqrt(diag(fisher_info))
      prop_sigma <- diag(prop_sigma)
    }
    else {
      prop_sigma <- 1 + pars/2
    }
  }
  I = diag(1e-8,length(pars))
  prop_sigma <- prop_sigma + I #makePositiveDefinite(prop_sigma)
  mu <- pars
  pi_X <- li_func(pars, full=T, sample=sample,...)
  k_X <- pars
  trace <- array(dim = c(iterations, length(pars)))
  llh = array(dim = iterations)
  accept = numeric(iterations)
  announce <- floor(seq(iterations/10, iterations, length.out = 10))
  eta = delta = list()
  for (i in 1:iterations) {
    k_Y <- mvnfast::rmvn(1, mu = k_X, sigma = prop_sigma)
    # cat('prop:',round(k_Y,2),' ')
    pi_Y <- li_func(k_Y, full=T, sample=sample,...)
    # cat('prop llh:',round(pi_Y,2),' ')
    a_X_Y = (pi_Y$ll) - (pi_X$ll)
    if (is.nan(a_X_Y))
      a_X_Y <- -Inf
    if (log(runif(1, 0, 1)) <= a_X_Y) {
      k_X = k_Y
      # cat('cur:',round(k_X,2),' ')
      pi_X = pi_Y
      accept[i] = 1
    } else{
      # if not accept recalc current state due to stochastic llh function
      # cat('cur:',round(k_X,2),' ')
      # pi_X = li_func(k_X, full=T,...)
      # resample w,v much faster than recalculating llh
      if(recalc){
        pi_X$eta = resample_w(pi_X$eta$w,flagp$Y.data$obs$trans,flagp$basis$obs$B,flagp$prior$ssq.prior.params)
        if(flagp$bias){
          pi_X$delta = resample_v(pi_X$delta,flagp$basis$obs$D,pi_X$eta$y.resid)
        }
        pi_X$ll = compute_ll(pi_X$theta,pi_X$ssq,pi_X$eta,pi_X$delta,sample,flagp,
                             flagp$prior$theta.prior,flagp$prior$theta.prior.params,flagp$prior$ssq.prior,flagp$prior$ssq.prior.params)
      }
      # if(sample.type=='logit')
      # pi_X$ll = pi_X$ll + logit.jacobian(t.curr) + log.jacobian(ssq.curr)
    }
    # cat('cur llh:',round(pi_X,2),'\n ')
    trace[i, ] <- k_X
    llh[i] = pi_X$ll
    eta[[i]] = pi_X$eta
    delta[[i]] = pi_X$delta

    if (i > adapt_par[1] && i%%adapt_par[2] == 0 && i < (adapt_par[4] * iterations)) {
      len <- floor(i * adapt_par[3]):i
      x <- trace[len, ]
      N <- length(len)
      p_sigma <- (N - 1) * var(x)/N ## SLOW
      # p_sigma <- makePositiveDefinite(p_sigma) ## can fail
      p_sigma = p_sigma + I
      if (!(0 %in% p_sigma))
        prop_sigma <- p_sigma
    }
    if (!quiet && i %in% announce)
      cat(paste0("  updating: ", i/iterations * 100, "%"," acceptance rate: ",round(mean(accept[1:i]),3),'\n'))
  }
  keep = seq(burn_in,iterations,thin)
  trace_all = trace
  trace <- trace[keep, ]
  llh = llh[keep]
  eta = eta[keep]
  delta = delta[keep]
  # if (length(pars) > 1) {
  #   theta_bar <- sapply(1:length(pars), function(x) {
  #     mean(trace[, x])
  #   })
  # }
  # else theta_bar <- mean(trace)
  accept_rate = mean(accept[burn_in:iterations])
  val <- list(t.samp = trace[,1:length(pars)-1,drop=F],ssq.samp=trace[,length(pars)], ll.samp = llh,
              eta = eta, delta = delta,
              prop.cov = prop_sigma, par_names = par_names,
              acpt.ratio = accept_rate, time=proc.time()[3]-ptm,n.samples=iterations,n.burn=burn_in,samp.all = trace_all)
  class(val) <- c("mcmc","list")
  return(val)
}

# MCMC sampling with theta on logit scale and s2 on log scale - in theory, sampling on the unconstrained space will improve things.
# in practice, covariance addaptation can be really sensitive in this space, blowing up so samples are always at
# the boundary
Metro_Hastings_Stochastic_Logit = function (flagp, li_func, pars, prop_sigma = NULL, par_names = NULL,
                                      iterations = 50000, burn_in = 1000, thin = 1, adapt_par = c(100, 20, 0.5, 0.75), recalc = TRUE, sample = TRUE, quiet = FALSE, ...)
{
  transform = function(pars,inv=F){
    p.t = length(pars)-1
    if(inv){
      pars[1:p.t] = invlogit(pars[1:p.t])
      pars[p.t+1] = exp(pars[p.t+1])
    } else{
      pars[1:p.t] = logit(pars[1:p.t])
      pars[p.t+1] = log(pars[p.t+1])
    }
    return(pars)
  }
  jacobian = function(pars,transformed=F){
    if(transformed)
      pars = transform(pars,inv=T) # pars need to be on regular scale
    sum(log(pars[1:p.t]*(1-pars[1:p.t]))) + log(pars[p.t+1])
  }
  ptm = proc.time()[3]
  if(is.null(sample))
    sample = as.logical(ifelse(flagp$bias,F,T))
  if (!is.finite(li_func(pars, sample=sample,...)))
    stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")
  if (is.null(par_names))
    par_names <- letters[1:length(pars)]
  if (!is.null(dim(prop_sigma))) {
    if ((dim(prop_sigma)[1] != length(pars) || dim(prop_sigma)[2] !=
         length(pars)) && !is.null(prop_sigma))
      stop("prop_sigma not of dimension length(pars) x length(pars)")
  }
  if (is.null(prop_sigma)) {
    if (length(pars) != 1) {
      fit <- optim(pars, li_func, control = list(fnscale = -1),
                   hessian = TRUE, opt=T, sample=sample,...)
      fisher_info <- solve(-fit$hessian)
      prop_sigma <- sqrt(diag(fisher_info))
      prop_sigma <- diag(prop_sigma)
    }
    else {
      prop_sigma <- 1 + pars/2
    }
  }
  I = diag(1e-8,length(pars))
  prop_sigma <- prop_sigma + I #makePositiveDefinite(prop_sigma)

  p.t = length(pars)-1
  pi_X <- li_func(pars, full=T, sample=sample,...)
  pi_X$ll = pi_X$ll + jacobian(pars,transformed=F)
  k_X <- transform(pars)
  trace <- array(dim = c(iterations, length(pars)))
  llh = array(dim = iterations)
  accept = numeric(iterations)
  announce <- floor(seq(iterations/10, iterations, length.out = 10))
  eta = delta = list()
  for (i in 1:iterations) {
    k_Y <- mvnfast::rmvn(1, mu = k_X, sigma = prop_sigma)
    pi_Y <- li_func(transform(k_Y,inv=T), full=T, sample=sample,...)
    pi_Y$ll = pi_Y$ll + jacobian(k_Y,transformed=T)

    a_X_Y = (pi_Y$ll) - (pi_X$ll)
    if (is.nan(a_X_Y))
      a_X_Y <- -Inf
    if (log(runif(1, 0, 1)) <= a_X_Y) {
      k_X = k_Y
      pi_X = pi_Y
      accept[i] = 1
    } else{
      # if not accept recalc current state due to stochastic llh function
      # resample w,v much faster than recalculating llh
      if(recalc){
        pi_X$eta = resample_w(pi_X$eta$w,flagp$Y.data$obs$trans,flagp$basis$obs$B,flagp$prior$ssq.prior.params)
        if(flagp$bias){
          pi_X$delta = resample_v(pi_X$delta,flagp$basis$obs$D,pi_X$eta$y.resid)
        }
        pi_X$ll = compute_ll(pi_X$theta,pi_X$ssq,pi_X$eta,pi_X$delta,sample,flagp,
                             flagp$prior$theta.prior,flagp$prior$theta.prior.params,
                             flagp$prior$ssq.prior,flagp$prior$ssq.prior.params) + jacobian(k_X,transformed = T)
      }
    }
    trace[i, ] <- k_X
    llh[i] = pi_X$ll
    eta[[i]] = pi_X$eta
    delta[[i]] = pi_X$delta

    if (i > adapt_par[1] && i%%adapt_par[2] == 0 && i < (adapt_par[4] * iterations)) {
      len <- floor(i * adapt_par[3]):i
      x <- trace[len, ]
      N <- length(len)
      p_sigma <- (N - 1) * var(x)/N ## SLOW
      # p_sigma <- makePositiveDefinite(p_sigma) ## can fail
      p_sigma = p_sigma + I
      if (!(0 %in% p_sigma))
        prop_sigma <- p_sigma
    }
    if (!quiet && i %in% announce)
      cat(paste0(". updating: ", i/iterations * 100, "%"," acceptance rate: ",round(mean(accept[1:i]),3),'\n'))
  }
  keep = seq(burn_in,iterations,thin)
  trace_all = t(apply(trace,1,transform,inv=T))
  trace <- trace_all[keep, ]
  llh = llh[keep]
  eta = eta[keep]
  delta = delta[keep]

  accept_rate = mean(accept[burn_in:iterations])
  val <- list(t.samp = trace[,1:length(pars)-1,drop=F],ssq.samp=trace[,length(pars)], ll.samp = llh,
              eta = eta, delta = delta,
              prop.cov = prop_sigma, par_names = par_names,
              acpt.ratio = accept_rate, time=proc.time()[3]-ptm,n.samples=iterations,n.burn=burn_in,samp.all = trace_all)
  class(val) <- c("mcmc","list")
  return(val)
}

#' @title FlaGP MCMC with joint proposal
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param flagp an \code{flagp} object.
#' @param t.init vector an inital value of the calibration parameters
#' @param prop.sd.init vector of initial proposal distribution sd's
#' @param n.samples integer number of mcmc samples to keep
#' @param n.burn integer number of initial samples to discard
#' @param adapt.par vector of parameters for addaptive metropolis. The first parameter it the first iteration to begin adaptation, the second is the frequency of adaptation, the third is the proportion of past samples to include for adaptation, and the fourth is the iteration to stop adaptation.
#' @param end.eta integer number of neighbors to use for laGP emulator
#' @param delta.method Model used for discrepancy. 'newGP' is a full GP using laGP default emperical bayes. 'rgasp' is a default RobustGaSP implementation and 'lagp' is an laGP predictor using the ALC criterion
#' @param start.delta initial neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param end.delta final neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param theta.prior prior for calibration parameters, 'beta' or 'unif'
#' @param theta.prior.params parameters for beta prior. Only applicable of theta.prior='beta'
#' @param ssq.prior.params gamma prior parameters for error variance facilitating conjugate updating
#' @param prev.samples Use in accordance with add_samples(). A matrix of samples can be passed in to be used for proposal covariance estimation.
#' @param recalc.llh Should the likelihood be recalculated at every iteration. Generally this should be False and only considered if sample=T. Recalculing the likelihood with a new sample from the emulator (and discrepancy) can be useful if mcmc mixing is a problem. This is only likely to happen if a discrepancy model is being used or if the emulator poorly fits the observations and no discrepancy model is used. Importantly, inference is approximate if recalculation is used, as posterior variance will be artificially increased.
#' @param sample Account for uncertainty in the emulator and discrepancy model by sampling rather than by including the variance in the likelihood. Default is sample=T when bias=F which works well because sampling variability in the emulator is usually small. For a biased calibration we recommend sample=F. This will increase the cost of likelihood calculation, but result in a much more stable mcmc.
#' @param verbose print status updates
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
mcmc = function(flagp,
                t.init=rep(.5,flagp$XT.data$p.t),
                ssq.init=.01,
                prop.cov=diag(c(rep(.005,flagp$XT.data$p.t),1e-6),flagp$XT.data$p.t+1),
                n.samples=10000,n.burn=1000,
                adapt.par = c(100,50,.5,n.burn),
                end.eta=50,
                delta.method='newGP',start.delta=6,end.delta=50,
                theta.prior='beta',theta.prior.params=c(2,2),
                ssq.prior='hcauchy',ssq.prior.params=c(.1),
                prev.samples=NULL,
                recalc.llh=F,
                sample=as.logical(ifelse(flagp$bias,F,T)),
                verbose=T)
{
  if(is.null(flagp$Y.data$n)){
    stop('Cannot perform calibration, no observed data in model.')
  }
  ptm = proc.time()
  # constants and initializations
  bias = flagp$bias
  p.t = flagp$XT.data$p.t

  t.curr = t.init
  ssq.curr = ssq.init

  llt = fit_model(t.curr,ssq.curr,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
                  theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
  llt$ll = llt$ll + logit.jacobian(t.curr) + log.jacobian(ssq.curr) # jacobian is added here because it may be hard to move from initial location otherwise (jacobian can be a large decrease in llh)

  llt.prev = llt
  ll.prop = list()

  # Storage
  t.store = matrix(nrow=n.samples,ncol=flagp$num$p.t)
  ssq.store = numeric(n.samples)
  eta.store = vector(mode='list',length=n.samples)
  if(bias){
    delta.store = vector(mode='list',length=n.samples)
  } else{
    delta.store = NULL
  }
  ll.store = accept = numeric(n.samples)

  # initialize adaptive parameters
  # s.d = 2.4^2/p.t
  eps = diag(sqrt(.Machine$double.eps),p.t+1)
  if(verbose)
    cat('MCMC Start #--',format(Sys.time(), "%a %b %d %X"),'--#\n')
  for(i in 1:n.samples){
    if(verbose){
      if(i %% 1000 == 0){
        cat('MCMC iteration',i,'#--',format(Sys.time(), "%a %b %d %X"),'--#\n')
      }
    }

    # propose new t
    prop = proposal(t.curr,ssq.curr,p.t,prop.cov+eps)
    t.prop = prop$t.prop
    ssq.prop = prop$ssq.prop

    ll.prop = fit_model(t.prop,ssq.prop,flagp,F,sample,end.eta,delta.method,start.delta,end.delta,F,
                        theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
    ll.prop$ll = ll.prop$ll + logit.jacobian(t.prop) + log.jacobian(ssq.prop)

    # we have to recompute this every time to account for the variability in the likelihood at t.curr,
    # otherwise, the chain does not mix well.
    if(recalc.llh & i>1){
      llt$eta = resample_w(llt$eta$w,flagp$Y.data$obs$trans,flagp$basis$obs$B,ssq.prior.params)
      if(flagp$bias){
        llt$delta = resample_v(llt$delta$v,llt$eta$y.resid,delta.method)
      }
      llt$ll = compute_ll(t.curr,ssq.curr,llt$eta,llt$delta,sample,flagp,theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)
      # if(sample.type=='logit')
      llt$ll = llt$ll + logit.jacobian(t.curr) + log.jacobian(ssq.curr)
      if(flagp$bias){
        llt$ssq.hat = llt$delta$ssq.hat
        llt$delta$ssq.hat = NULL
      } else{
        llt$ssq.hat = llt$eta$ssq.hat
        llt$eta$ssq.hat = NULL
      }
    }

    # 4. Accept or reject t
    accept.prob =  exp(ll.prop$ll - llt$ll)
    if(runif(1) < accept.prob){
      t.curr = t.prop
      ssq.curr = ssq.prop
      if(!recalc.llh){
        # accepted, but not recalculating, so update old llt to ll.prop
        llt = ll.prop
      } else{
        # accepted and recalculating, update llt.prev
        llt.prev = ll.prop
      }
      # ssq.store[i] = ll.prop$ssq.hat
      ssq.store[i] = ssq.prop
      eta.store[[i]] = ll.prop$eta
      if(bias){delta.store[[i]] = ll.prop$delta}
      ll.store[i] = ll.prop$ll
      accept[i] = 1
    } else{
      if(recalc.llh){
        llt = llt.prev
      }
      # ssq.store[i] = llt$ssq.hat
      eta.store[[i]] = llt$eta
      if(bias){
        delta.store[[i]] = llt$delta
      }
      ll.store[i] = llt$ll
    }
    t.store[i,] = t.curr
    ssq.store[i] = ssq.curr

    # adapt proposal distribution
    if(i>=adapt.par[1] & i<=adapt.par[4] & i %% adapt.par[2] == 0){
      interval = max(1,i-floor(adapt.par[3]*i)):i
      t.interval = t.store[interval,,drop=F]
      ssq.interval = ssq.store[interval]
      lambda = ifelse(i==adapt.par[1],1,update$lambda)
      update = update_tuning_mv(i, mean(accept[interval]), lambda=lambda, cbind(log(t.interval/(1-t.interval)),ssq.interval),
                                prop.cov/lambda)
      prop.cov = update$lambda * update$Sigma_tune
    }
  }
  mcmc.time = proc.time() - ptm

  return.ids = (n.burn+1):n.samples
  returns = list(t.samp     = t.store[return.ids,,drop=F],
                 ssq.samp   = ssq.store[return.ids],
                 ll.samp    = ll.store[return.ids],
                 eta        = eta.store[return.ids],
                 delta      = delta.store[return.ids],
                 prop.cov   = prop.cov,
                 acpt.ratio = mean(accept[return.ids]),
                 time       = mcmc.time,
                 n.samples  = n.samples,
                 n.burn     = n.burn)
  class(returns) <- c("mcmc", class(returns))
  return(returns)
}
logit.jacobian = function(theta){
  sum(log(theta)+log(1-theta))
}
log.jacobian = function(ssq){
  log(ssq)
}

proposal = function(t.curr,ssq.curr,p.t,prop.cov){
  # propose theta on logit scale and error variance on log scale
  par.prop = mvnfast::rmvn(1,c(log(t.curr/(1-t.curr)),log(ssq.curr)), prop.cov)
  t.prop = exp(par.prop[1:p.t])/(1+exp(par.prop[1:p.t]))

  # NA's can happen when proposal covariance gets too big exp(huge)=Inf
  t.prop[is.na(t.prop)] = 1
  t.prop[t.prop<.Machine$double.eps] = .Machine$double.eps

  return(list(t.prop=t.prop,
              ssq.prop=exp(par.prop[p.t+1])))
}
# Function to add samples to existing mcmc object, starting from the last interation.
# mcmc() has been updated, this function needs to be as well to reflect that
add_samples = function(flagp,mcmc,n.samples=100,
                       end.eta=50,sample=F,recalc.llh=F,
                       delta.method='newGP',start.delta=6,end.delta=50,
                       adapt.par = c(100,50,50,1000),
                       theta.prior='beta',theta.prior.params=c(2,2),
                       ssq.prior='hcauchy',ssq.prior.params=c(.5),
                       verbose=T){

  add.samples = mcmc(flagp,
                     t.init=mcmc$t.samp[nrow(mcmc$t.samp),],
                     ssq.init=mcmc$ssq.samp[length(mcmc$ssq.samp)],
                     prop.cov=mcmc$prop.cov,
                     n.samples=n.samples,n.burn=0,adapt.par=adapt.par,
                     end.eta=end.eta,delta.method=delta.method,start.delta=start.delta,end.delta=end.delta,
                     theta.prior=theta.prior,
                     ssq.prior=ssq.prior,ssq.prior.params=ssq.prior.params,
                     prev.samples=mcmc$t.samp,recalc.llh=recalc.llh,
                     sample=sample,verbose=verbose)

  mcmc$t.samp = rbind(mcmc$t.samp,add.samples$t.samp)
  mcmc$ssq.samp = c(mcmc$ssq.samp,add.samples$ssq.samp)
  mcmc$ll.samp = c(mcmc$ll.samp,add.samples$ll.samp)
  mcmc$eta = append(mcmc$eta,add.samples$eta)
  mcmc$delta = append(mcmc$delta,add.samples$delta)

  # acceptance ratio is weighted avg of number of samples
  prev.samples = mcmc$n.samples-mcmc$n.burn
  mcmc$acpt.ratio = (prev.samples*mcmc$acpt.ratio + n.samples*add.samples$acpt.ratio)/(prev.samples+n.samples)
  mcmc$time = mcmc$time + add.samples$time
  mcmc$n.samples = mcmc$n.samples + n.samples
  return(mcmc)
}

#' @title Parallel wrapper for mcmc function
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param flagp an \code{flagp} object.
#' @param t.init vector an inital value of the calibration parameters
#' @param prop.sd.init matrix initial covariance matrix for the proposal distribution
#' @param n.samples integer number of mcmc samples to keep
#' @param n.burn integer number of initial samples to discard
#' @param update.every integer sample frequency that proposal covariance should be updated, only working for update.every=1 currently
#' @param stop.update integer sample index where proposal covariance updates should stop
#' @param end.eta integer number of neighbors to use for laGP emulator
#' @param delta.method Model used for discrepancy. 'newGP' is a full GP using laGP default emperical bayes. 'rgasp' is a default RobustGaSP implementation and 'lagp' is an laGP predictor using the ALC criterion
#' @param start.delta initial neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param end.delta final neighborhood size for bias model. Only applicable if delta.method='lagp'
#' @param theta.prior prior for calibration parameters, 'beta' or 'unif'
#' @param ssq.prior.params gamma prior parameters for error variance facilitating conjugate updating
#' @param prev.samples Use in accordance with add_samples(). A matrix of samples can be passed in to be used for proposal covariance estimation.
#' @param verbose print status updates
#' @param n.chains number of parallel chains to run
#' @param n.cores number of cores for parallel processing
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
parallel_mcmc = function(flagp,t.init=rep(.5,flagp$XT.data$p.t),
                         prop.cov=diag((.5/3)^2,flagp$XT.data$p.t),
                          n.samples=10000,n.burn=1000,
                          update.every=50,stop.update=1000,
                          end.eta=50,
                          delta.method='newGP',start.delta=6,end.delta=50,
                          theta.prior='beta',theta.prior.params=c(2,2),
                          ssq.prior.params=c(1,1e-3),prev.samples=NULL,
                          recalc.llh=F,
                          verbose=T,n.chains=2,n.cores=2)
{
  ptm = proc.time()

  cl = parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)

  # constants and initializations
  bias = flagp$bias
  p.t = flagp$XT.data$p.t
  if(is.null(t.init) | nrow(t.init)!=n_chains){
    t.init = tgp::lhs(n.chains,matrix(rep(c(0,1),n.chains),nrow=n.chains,byrow=T))
  }

  # do parallel mcmc chains
  results = foreach::foreach(i=1:n.chains) %dopar% mcmc(flagp, t.init[i,], prop.cov, n.samples, n.burn, update.every, stop.update,
                                                        end.eta, delta.method, start.delta, end.delta,
                                                        theta.prior, theta.prior.params,
                                                        ssq.prior.params, prev.samples, recalc.llh, verbose)

  parallel::stopCluster(cl)
  # combine results
  results$t.samp = c()
  results$ssq.samp = c()
  for(i in 1:length(results)){
    results$t.samp = rbind(results$t.samp,results[[i]]$t.samp)
    results$ssq.samp = c(results$ssq.samp,results[[i]]$ssq.samp)
  }
  class(results) <- c("mcmc", class(results))
  return(results)
}

#' this function updates a block Gaussian random walk proposal Metropolis-Hastings tuning parameters
#' @param k is a positive integer that is the current MCMC iteration.
#' @param accept is a positive number between 0 and 1 that represents the batch accpetance rate. The target of this adaptive tuning algorithm is an acceptance rate of between 0.44 (a univariate update) and 0.234 (5+ dimensional upate).
#' @param lambda is a positive scalar that scales the covariance matrix Sigma_tune and diminishes towards 0 as k increases.
#' @param batch_samples is a \eqn{batch.size \times d}{batch.size x d} dimensional matrix that consists of the batch samples of the d-dimensional parameter being sampled.
#' @param Sigma_tune is a \eqn{d \times d}{d x d} positive definite covariance matrix of the batch samples used to generate the multivariate Gaussian random walk proposal.
#'
#' @keywords internal
update_tuning_mv <- function(k, accept, lambda, batch_samples,
                             Sigma_tune) {
  # target acceptance rate by dimension of parameter vector, min of .234
  arr <- c(0.44, 0.35, 0.32, 0.25, 0.234)
  d = ncol(batch_samples)
  dimension <- min(d,5)
  batch_size <- nrow(batch_samples)
  optimal_accept <- arr[dimension]
  times_adapted <- floor(k / batch_size)
  gamma1 <- 1.0 / ((times_adapted + 3.0)^0.8)
  gamma2 <- 10.0 * gamma1
  adapt_factor <- exp(gamma2 * (accept - optimal_accept))
  ## update the MV scaling parameter
  lambda_out <- lambda * adapt_factor
  ## center the batch of MCMC samples
  batch_samples_tmp <- batch_samples
  for (j in 1:d) {
    mean_batch = mean(batch_samples[, j])
    for (i in 1:batch_size) {
      batch_samples_tmp[i, j] <- batch_samples[i, j] - mean_batch
    }
  }
  Sigma_tune_out <- Sigma_tune + gamma1 *
    (t(batch_samples_tmp) %*% batch_samples_tmp / (batch_size-1.0) - Sigma_tune)
  accept_out <- 0.0
  batch_samples_out <- matrix(0, batch_size, d)
  return(
    list(
      accept          = accept_out,
      lambda          = lambda_out,
      batch_samples   = batch_samples_out,
      Sigma_tune      = Sigma_tune_out
    )
  )
}

Metro_Hastings_Stochastic_Logit_OneAtATime = function (
    flagp,
    li_func,
    pars,
    prop_sigma = NULL,
    par_names = NULL,
    iterations = 50000,
    burn_in = 1000,
    thin = 1,
    adapt_par = c(100, 20, 0.5, 0.75),
    recalc = TRUE,
    sample = TRUE,
    quiet = FALSE,
    ...
) {
  ##
  ## 1) Define transformations and Jacobian
  ##
  transform = function(pars, inv = FALSE){
    p.t = length(pars) - 1
    if(inv) {
      # back-transform
      pars[1:p.t] = invlogit(pars[1:p.t])
      pars[p.t+1] = exp(pars[p.t+1])
    } else {
      # forward transform
      pars[1:p.t] = logit(pars[1:p.t])
      pars[p.t+1] = log(pars[p.t+1])
    }
    return(pars)
  }

  jacobian = function(pars, transformed = FALSE){
    # Returns sum of log of absolute value of the derivative
    # for each parameter dimension
    p.t = length(pars) - 1
    if(transformed) {
      # 'pars' is currently on *transformed* scale
      # so we invert it to original scale before computing
      pars = transform(pars, inv = TRUE)
    }
    # Now 'pars' are on the original scale:
    # partial derivatives:
    # d/d(logit(p)) = 1/[p(1-p)],   so log(|d/d(logit(p))|) = -log([p(1-p)])
    # d/d(log(sigma)) = 1/sigma,    so log(|d/d(log(sigma))|) = -log(sigma)
    # => But we often incorporate these as + log(Jacobian) or - log(Jacobian)
    # depending on how you structured li_func.  In your original code, you do:
    #   sum(log(p_i*(1-p_i))) + log(sigma)
    # so let's keep that approach.
    jac = sum(log(pars[1:p.t] * (1 - pars[1:p.t]))) + log(pars[p.t+1])
    return(jac)
  }

  ##
  ## 2) Initial checks and setup
  ##
  ptm = proc.time()[3]
  if(is.null(sample))
    sample = as.logical(ifelse(flagp$bias, TRUE, FALSE))

  # Evaluate log-likelihood at starting values
  if (!is.finite(li_func(pars, sample = sample, ...))) {
    stop("Seed parameter values <pars> are not in the defined parameter space.
         Try new starting values for <pars>.")
  }

  if (is.null(par_names))
    par_names <- letters[1:length(pars)]

  # If user has provided a proposal sigma, check dims
  if (!is.null(dim(prop_sigma))) {
    if ((dim(prop_sigma)[1] != length(pars) || dim(prop_sigma)[2] != length(pars))) {
      stop("prop_sigma not of dimension length(pars) x length(pars)")
    }
  }

  ## 3) If no proposal covariance is given, do a quick optimization-based guess
  if (is.null(prop_sigma)) {
    if (length(pars) > 1) {
      fit <- optim(pars, li_func, control = list(fnscale = -1),
                   hessian = TRUE, opt = TRUE, sample=sample, ...)
      fisher_info <- solve(-fit$hessian)
      stdevs <- sqrt(diag(fisher_info))
      prop_sigma <- diag(stdevs)
    } else {
      # Single-parameter guess
      prop_sigma <- matrix(1 + pars/2, nrow = 1)
    }
  }
  # Force positive-definiteness by adding small identity
  I = diag(1e-8, length(pars))
  prop_sigma <- prop_sigma + I

  p.t = length(pars) - 1

  # Evaluate 'pi_X' at the starting point, on original scale:
  pi_X <- li_func(pars, full = TRUE, sample = sample, ...)
  # Add Jacobian from the untransformed 'pars':
  pi_X$ll = pi_X$ll + jacobian(pars, transformed = FALSE)

  # 'k_X' is in the transformed space
  k_X <- transform(pars)

  # Storage
  trace <- matrix(NA, nrow = iterations, ncol = length(pars))
  llh   <- numeric(iterations)
  accept = numeric(iterations)  # acceptance count per *iteration* (summing over params)

  # Possibly store any states from pi_X
  eta   <- vector("list", iterations)
  delta <- vector("list", iterations)

  # For printing progress
  announce <- floor(seq(iterations/10, iterations, length.out = 10))

  ##
  ## 4) Main iteration loop
  ##
  for (iter in seq_len(iterations)) {

    # For each parameter dimension, do a one-at-a-time update
    n_accepted = 0  # track how many accepted in this iteration

    for (j in seq_along(k_X)) {

      # Propose new value for dimension j
      k_Y = k_X
      # sample from Normal(k_X[j], sqrt(prop_sigma[j,j]))
      # you can use rnorm or anything else
      k_Y[j] = rnorm(1, mean = k_X[j], sd = sqrt(prop_sigma[j,j]))

      # Evaluate new log-likelihood
      pi_Y <- li_func(transform(k_Y, inv = TRUE), full = TRUE, sample = sample, ...)
      # Add Jacobian on transformed scale
      pi_Y$ll = pi_Y$ll + jacobian(k_Y, transformed = TRUE)

      # MH acceptance ratio (log scale)
      a_X_Y = pi_Y$ll - pi_X$ll

      if (is.nan(a_X_Y)) a_X_Y <- -Inf

      # Accept/reject
      if (log(runif(1)) <= a_X_Y) {
        # Accept
        k_X  = k_Y
        pi_X = pi_Y
        n_accepted = n_accepted + 1
      } else {
        # Reject
        # If recalc = TRUE, we re-sample some latent variables in pi_X
        if (recalc) {
          pi_X$eta = resample_w(pi_X$eta$w,
                                flagp$Y.data$obs$trans,
                                flagp$basis$obs$B,
                                flagp$prior$ssq.prior.params)
          if (flagp$bias) {
            pi_X$delta = resample_v(pi_X$delta,
                                    flagp$basis$obs$D,
                                    pi_X$eta$y.resid)
          }
          pi_X$ll = compute_ll(pi_X$theta, pi_X$ssq, pi_X$eta, pi_X$delta,
                               sample, flagp,
                               flagp$prior$theta.prior,
                               flagp$prior$theta.prior.params,
                               flagp$prior$ssq.prior,
                               flagp$prior$ssq.prior.params)
          # Add jacobian for the *current* k_X
          pi_X$ll = pi_X$ll + jacobian(k_X, transformed = TRUE)
        }
      }
    }  # end loop over parameters

    # Store
    trace[iter, ] = k_X
    llh[iter]     = pi_X$ll
    eta[[iter]]   = pi_X$eta
    delta[[iter]] = pi_X$delta

    accept[iter]  = n_accepted / length(k_X)

    ##
    ## Adaptation step
    ##
    if (iter > adapt_par[1] &&
        (iter %% adapt_par[2]) == 0 &&
        iter < (adapt_par[4] * iterations)) {

      len <- floor(iter * adapt_par[3]):iter
      x   <- trace[len, ]
      N   <- length(len)

      # Compute diagonal from var(x) on the transformed scale
      # var(x) returns a full covariance; we take diag and rebuild
      # p_var = diag(var(x) * (N-1)/N)  # unbiased estimate
      # p_var = p_var + I   # ensure positivity
      p_var = apply(x,2,var)
      prop_sigma = diag(p_var)
    }

    # Progress printing
    if(!quiet && iter %in% announce) {
      cat(paste0(". iteration: ", iter,
                 " (", round(iter/iterations * 100, 2), "%), ",
                 " acceptance rate: ",
                 round(mean(accept[1:iter]), 3), "\n"))
    }
  }

  ##
  ## 5) Post-processing
  ##
  keep       = seq(burn_in, iterations, by = thin)
  trace_all  = t(apply(trace, 1, transform, inv = TRUE))
  trace_keep = trace_all[keep, , drop = FALSE]
  llh_keep   = llh[keep]
  eta_keep   = eta[keep]
  delta_keep = delta[keep]

  accept_rate = mean(accept[burn_in:iterations])

  val <- list(
    t.samp     = trace_keep[, 1:(length(pars)-1), drop=FALSE],
    ssq.samp   = trace_keep[, length(pars)],
    ll.samp    = llh_keep,
    eta        = eta_keep,
    delta      = delta_keep,
    prop.cov   = prop_sigma,
    par_names  = par_names,
    acpt.ratio = accept_rate,
    time       = proc.time()[3] - ptm,
    n.samples  = iterations,
    n.burn     = burn_in,
    samp.all   = trace_all
  )
  class(val) <- c("mcmc","list")
  return(val)
}
