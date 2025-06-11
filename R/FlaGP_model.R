fit_model = function(theta,ssq,flagp,
                     lite=T,sample=T,
                     end.eta=50,
                     delta.method='newGP',start.delta=6,end.delta=50,
                     negll=F,
                     theta.prior='beta',theta.prior.params=c(2,2),
                     ssq.prior='hcauchy',ssq.prior.params=c(.5),map=F){
  # Predict from emulator at [X.obs,theta]
  eta = FlaGP:::fit_eta(theta,flagp,sample,end.eta,ssq.prior.params,map=map)

  if(flagp$bias){
    y.resid = eta$y.resid
    D = flagp$basis$obs$D
    XT = flagp$XT.data
    delta = FlaGP:::fit_delta(y.resid,D,XT,sample=sample,delta.method=delta.method,start=start.delta,end=end.delta,ssq.prior.params = ssq.prior.params,map=map)
  } else{
    delta = NULL
  }

  ll = FlaGP:::compute_ll(theta,ssq,eta,delta,sample,flagp,theta.prior,theta.prior.params,ssq.prior,ssq.prior.params)

  if(lite){
    return(ifelse(negll,-ll,ll))
  } else{
    # if(flagp$bias){
    #   ssq.hat = delta$ssq.hat
    #   eta$ssq.hat = NULL; delta$ssq.hat = NULL
    # } else{
    #   ssq.hat = eta$ssq.hat
    #   eta$ssq.hat = NULL
    # }

    return(list(ll = ifelse(negll,-ll,ll),
                theta = theta,
                ssq = ssq,
                eta = eta,
                delta = delta))
  }
}
# wrapper function for fit model that has theta and ssq contained in the first argument for use with optimization function
fit_model_map = function(param,flagp,
                         end.eta=50,
                         delta.method='newGP',start.delta=6,end.delta=50,
                         theta.prior='beta',theta.prior.params=c(2,2),
                         ssq.prior='hcauchy',ssq.prior.params=c(.5),loglogit=F){

  if(loglogit){
    theta = LaplacesDemon::invlogit(param[1:flagp$XT.data$p.t]) # theta passed on logit scale for optimization
    ssq = exp(param[flagp$XT.data$p.t+1]) # passed in on log scale
  } else{
    theta = param[1:flagp$XT.data$p.t]
    ssq = param[flagp$XT.data$p.t+1]
  }

  if(!all(theta < 1 & theta > 0) | ssq<=0){
    return(Inf) # minimization problem
  }
  FlaGP:::fit_model(theta=theta,ssq=ssq,flagp=flagp,
                    lite=T,sample=F,
                    end.eta=end.eta,
                    delta.method=delta.method,start.delta=start.delta,end.delta=end.delta,
                    negll=T,theta.prior=theta.prior,theta.prior.params=theta.prior.params,
                    ssq.prior.params=ssq.prior.params,map=T)
}
# fit modular emulator with t=theta
fit_eta = function(theta,flagp,sample=T,end=50,ssq.prior.params,map)
{

  n = flagp$Y.data$n
  n.pc = flagp$basis$sim$n.pc
  p.t = flagp$XT.data$p.t
  p.x = flagp$XT.data$p.x

  # theta was sampled on 0-1, so we need to stretch and compress it by our lengthsale estimates
  theta.sc = lapply(1:n.pc, function(jj) theta/sqrt(flagp$lengthscales$T[[jj]]))
  theta.sc.rep = lapply(1:n.pc, function(jj) matrix(theta.sc[[jj]],nrow=n,ncol=p.t,byrow = T))

  # make prediction matrix from X.obs and theta
  if(p.t>1){
    if(p.x>0){
      XT.sc = lapply(1:n.pc, function(jj) cbind(flagp$SC.inputs$X.obs[[jj]],theta.sc.rep[[jj]]))
    } else{
      XT.sc = theta.sc.rep
    }
  } else{
    XT.sc = lapply(1:n.pc, function(jj) cbind(flagp$SC.inputs$X.obs[[jj]],theta.sc.rep[[jj]]))
  }

  w = aGPsep_SC_mv(X=flagp$SC.inputs$XT.sim,
                   Z=flagp$basis$sim$V.t,
                   XX=XT.sc,
                   end=end,
                   sample=sample,g=flagp$lengthscales$g)

  # residuals
  y.pred = w_to_y(w$sample,flagp$basis$obs$B)
  y.resid = flagp$Y.data$obs$trans - y.pred

  # conjugate posterior
  # rtr = apply(y.resid,2,function(x) t(x) %*% x)
  # n.s2 = prod(dim(flagp$Y.data$obs$trans))
  # if(!map){
  #   ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)
  # } else{
  #   # mean of inv gamma for MAP
  #   ssq.hat = (ssq.prior.params[2]+sum(rtr)/2)/(ssq.prior.params[1]+n.s2/2-1)
  # }

  return(list(w=w,y.resid=y.resid))#,ssq.hat=ssq.hat)))
}
# In mcmc we has previously recalculated the likelihood at the current theta to improve mixing.
# This is cumbersome when all we really needed was to resample the emulator and discrepancy model
resample_w = function(w,y.obs.trans,B.obs,ssq.prior.params){
  w$sample = rt(prod(dim(w$mean)),w$df) * sqrt(w$scale) + w$mean
  y.pred = w_to_y(w$sample,B.obs)
  y.resid = y.obs.trans - y.pred

  # conjugate posterior
  # rtr = apply(y.resid,2,function(x) t(x) %*% x)
  # n.s2 = prod(dim(y.obs.trans))
  # ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)

  return(list(w=w,y.resid=y.resid))#,ssq.hat=ssq.hat))
}

# fit delta model to residuals using basis vectors in D
fit_delta = function(y.resid,D,XT.data,sample,delta.method,start,end,ssq.prior.params,map){
  returns = list()
  v = NULL
  n.pc = ncol(D)
  # Center y.resid first, because we want v to be mean zero when far away from training data
  #y.resid.mean = apply(y.resid,1,mean)
  #y.resid = y.resid - y.resid.mean # center y.resid - zero mean GP
  #y.resid.sd = sd(y.resid)
  #y.resid = y.resid / y.resid.sd
  V.t = get_basis(y.resid,B=D)$V.t
  v$df = ncol(V.t)
  #V.t.mean = apply(V.t,2,mean)
  #V.t = sweep(V.t,2,V.t.mean)
  #V.t.sd = sd(V.t)
  #V.t = V.t / V.t.sd
  if(delta.method=='lagp'){ # use laGP for a faster bias model
    v = append(v,aGPsep_SC_mv(X=XT.data$obs$X$trans,
                              Z=V.t,
                              XX=XT.data$obs$X$trans,
                              start=start,
                              end=end,
                              bias=T))
    v$df = end
  } else if(delta.method=='newGP'){       # use full GP for the bias model
    GPs = vector(mode='list',length=n.pc)
    mle = vector(mode='list',length=n.pc)
    for(k in 1:n.pc){
      if(ncol(y.resid)>=5){
        # darg fails for one observation and generally these defaults don't seem to work well for small data problems
        d = laGP::darg(d=list(mle=T),X=XT.data$obs$X$trans)
        g = garg(g=list(mle=T),y=V.t[k,])
      } else{
        # we don't want too wiggly of a discrepancy so d should be large
        d = list(mle=T,start=5,min=1,max=10,ab=c(5,1))
        # we want to encourage it to be near zero, so what to do with g?
        g = list(mle=T,start=.01,min=1e-8,max=1,ab=c(0,0))
      }

      # Our responses have changed w.r.t to inputs, so we can't not use d=1 anymore
      GPs[[k]] <- laGP::newGPsep(X=XT.data$obs$X$trans,
                             Z=V.t[k,],
                             d=d$start, g=g$start, dK=TRUE)
      cmle <- laGP::jmleGPsep(GPs[[k]],
                        drange=c(d$min, 10*d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      mle[[k]] = cmle
      pred = laGP::predGPsep(GPs[[k]], XX = XT.data$obs$X$trans, lite=T)

      v$mean = rbind(v$mean, pred$mean)
      pred$s2[pred$s2<0]=0
      v$var = rbind(v$var,pred$s2)
      laGP::deleteGPsep(GPs[[k]])
    }
    v$scale = v$var*(v$df-2)/v$df # convert variance to scale of student t
    returns$mle=mle
  } else if(delta.method=='rgasp'){
    returns$model = vector(mode='list',length=n.pc)
    for(k in 1:n.pc){
      invisible(capture.output(returns$model[[k]] <- RobustGaSP::rgasp(XT.data$obs$X$trans, V.t[k,],nugget=0,nugget.est=T,num_initial_values = 1)))
      pred = RobustGaSP::predict(returns$model[[k]],XT.data$obs$X$trans)

      v$mean = rbind(v$mean, pred$mean)
      v$var = rbind(v$var, pred$sd^2)
    }
  } else if(delta.method=='mlegp'){
    returns$model = vector(mode='list',length=n.pc)
    for(k in 1:n.pc){
      invisible(capture.output(returns$model[[k]] <- mlegp::mlegp(XT.data$obs$X$trans, V.t[k,],simplex.ntries = 1,verbose = 0, min.nugget = sqrt(.Machine$double.eps))))
      pred = predict(returns$model[[k]], se.fit = T)

      v$mean = rbind(v$mean, t(pred$fit))
      v$var = rbind(v$var, t(pred$se.fit)^2)
    }
  } else if(delta.method=='homgp'){
    returns$model = vector(mode='list',length=n.pc)
    for(k in 1:n.pc){
      returns$model[[k]] <- hetGP::mleHomGP(XT.data$obs$X$trans, V.t[k,])
      pred = predict(object = returns$model[[k]], x=XT.data$obs$X$trans)

      v$mean = rbind(v$mean, pred$mean)
      v$var = rbind(v$var, pred$sd2)
    }
  } else{
    stop('delta.method must be either lagp, newGP, rgasp, mlegp, or homgp.')
  }
  # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
  if(sample){
    # sampling via cholesky where chol = sqrt(var)
    if(delta.method %in% c('lagp','newGP')){
      if(all(v$scale>0)){
        v$sample = (rt(prod(dim(v$mean)),v$df) * sqrt(v$scale)) + v$mean
      } else{
        v$sample = (rnorm(prod(dim(v$mean))) * sqrt(v$var)) + v$mean
      }
    } else{
      v$sample = (rnorm(prod(dim(v$mean))) * sqrt(v$var)) + v$mean
    }
  } else{
    v$sample = v$mean
  }
  returns$v = v
  returns$V.t = V.t
  returns$y.resid = y.resid - w_to_y(v$sample,D)
  #returns$y.resid = (y.resid + y.resid.mean) - (w_to_y(v$sample,D) + y.resid.mean)
  #returns$v$y.resid.mean = y.resid.mean

  # rtr = apply(returns$y.resid,2,function(x) t(x) %*% x)
  # n.s2 = prod(dim(returns$y.resid))
  # if(!map){
  #   returns$ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)
  # } else{
  #   # mean of inv gamma for MAP
  #   returns$ssq.hat = (ssq.prior.params[2]+sum(rtr)/2)/(ssq.prior.params[1]+n.s2/2-1)
  # }
  returns$method = delta.method
  return(returns)
}

resample_v = function(delta,D,eta.y.resid){##ssq.prior.params,delta.method){
  # resample v
  if(delta$method %in% c('lagp','newGP')){
    if(all(delta$v$scale>0)){
      delta$v$sample = (rt(prod(dim(delta$v$mean)),delta$v$df) * sqrt(delta$v$scale)) + delta$v$mean
    } else{
      delta$v$sample = (rnorm(prod(dim(delta$v$mean))) * sqrt(delta$v$var)) + delta$v$mean
    }
  } else{
    delta$v$sample = (rnorm(prod(dim(delta$v$mean))) * sqrt(delta$v$var)) + delta$v$mean
  }
  # recompute residuals
  delta$y.resid = eta.y.resid - w_to_y(delta$v$sample,D)# + v$y.resid.mean)

  # conjugate posterior
  # rtr = apply(y.resid,2,function(x) t(x) %*% x)
  # n.s2 = prod(dim(y.resid))
  # ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)

  return(delta)#,ssq.hat=ssq.hat))
}

# FUNCTION: multivariate aGPsep for use with stretched and compressed inputs only
{## Parameters:
# required
# X : list of length n.pc which contains stretched and compressed training input matrices
# Z : list of length n.pc which contains training response vectors
# XX: list of length n.pc with contains stretched and compressed prediction input matrices
# optional
# start: integer, initial size of laPG neighborhood
# end: integer, final size of laGP neighborhood
# g: double, laGP nugget
## Returns:
# lagp_fit: list of n.pc outputs from aGPsep
# mean: matrix (n.pc x n) of prediction means
# scale: matrix (n.pc x n) of prediction scale parameters (student t with end degrees of freedom)
## Calls:
# aGPsep from laGP package
# w_to_y to transform w -> y
## Called by:
# mv_calib.bias: biased calibration
# mv_calib.nobias: unbiased calibration
# mv.em.pred: laGP emulator prediction function at new inputs
}
aGPsep_SC_mv = function(X, Z, XX, g, start=6, end=50, bias=F, sample=F, predvar=T, parallel = F){
  select_maximin_subset_indices <- function(points, subset_size) {
    n_points <- nrow(points)

    if(subset_size > n_points){
      stop("Subset size cannot exceed the number of available points.")
    }

    # Initialize with the first point
    selected_indices <- c(1)
    remaining_indices <- setdiff(1:n_points, selected_indices)

    for(i in 2:subset_size){
      # Compute the minimum distance of each remaining point to the selected points
      min_distances <- sapply(remaining_indices, function(idx){
        min(sqrt(rowSums((t(points[selected_indices, ]) - points[idx, ])^2)))
      })

      # Select the point with the maximum of these minimum distances
      next_idx <- remaining_indices[which.max(min_distances)]
      selected_indices <- c(selected_indices, next_idx)
      remaining_indices <- setdiff(remaining_indices, next_idx)
    }

    return(selected_indices)
  }
  n.pc = nrow(Z)
  n.XX = ifelse(bias,nrow(XX),nrow(XX[[1]]))
  mean = scale = krigvar = array(dim=c(n.pc,n.XX))
  if(bias){
    # if we use laGP for the bias, we cannot use stretched an compressed inputs anymore so default to alc criterion
    lagp_fit = lapply(1:n.pc,function(i) laGP::aGPsep(X = X,
                                                Z = Z[i,],
                                                XX = XX,
                                                start = start,
                                                end = min(end,ncol(Z)-1),
                                                method = 'alc',
                                                verb=0))
    # TODO: I don't think this is working, we need to take the results from lagp_fit and make a mean/scale matrix
  } else{
    # SUPER SLOW
    # lagp_pred = lapply(1:n.pc, function(i) laGP::aGPsep(X[[i]],Z[i,],XX[[i]],method = 'nn',d=list(start=1,mle=F),g=g,verb = F,end = end))
    # mean = t(sapply(1:n.pc, function(i) lagp_pred[[i]]$mean))
    # var = t(sapply(1:n.pc, function(i) lagp_pred[[i]]$var))

    if(end<nrow(X[[1]])){
      nn.indx = lapply(1:n.pc, function(i) FNN::get.knnx(X[[i]],XX[[i]],end)$nn.index)
    } else{
      nn.indx = lapply(1:n.pc,function(i) matrix(rep(seq(1:nrow(X[[i]])),n.XX),nrow=n.XX,byrow = T))
    }

    if(!predvar){
      # If we don't need the prediction variance, we can shortcut the mean by computing the Cholesky of [X,XX] together
      for(j in 1:n.XX){
        for(i in 1:n.pc){
          GP = GP_fit_isotropic(rbind(X[[i]][nn.indx[[i]][j,],],XX[[i]][j,]),
                                d=1,g[i],lite=T)
          tKchol = t(GP$Kchol)
          # cross-cov: K[end+1,1:end]
          # train-cov: K[1:end,1:end]
          # traindata: Z[i,nn.indx[[i]][j,]]
          # Kt %*% KiY
          mean[i,j] = tKchol[end+1,1:end]%*%forwardsolve(tKchol[1:end,1:end],Z[i,nn.indx[[i]][j,]])
        }
      }
    } else{
      if(parallel){
        Xadd = NULL
        parallel_pred = function(i){
          MVK = array(dim=c(1,n.XX,3))
          for(j in 1:n.XX){
            lagp_fit = GP_fit_isotropic(rbind(X[[i]][nn.indx[[i]][j,],,drop=F],X[[i]][Xadd[[i]],]),
                             d=1,g[i],
                             c(Z[i,nn.indx[[i]][j,]],Z[i,Xadd[[i]]]),
                             lite=F)
            lagp_pred = GP_predict(lagp_fit,
                                   rbind(X[[i]][nn.indx[[i]][j,],,drop=F],X[[i]][Xadd[[i]],]),
                                   XX[[i]][j,,drop=F],
                                   c(Z[i,nn.indx[[i]][j,]],Z[i,Xadd[[i]]]),
                                   predvar=T)
            MVK[1,j,1] = lagp_pred$mean
            MVK[1,j,2] = lagp_pred$s2
            MVK[1,j,3] = lagp_pred$krigvar
          }
          return(MVK)
        }
        combine_matrices <- function(x, y) {
          abind::abind(x, y, along = 1)  # Stack along the third dimension
        }
        MVK = foreach::foreach(i=1:n.pc,.combine = combine_matrices) %dopar% parallel_pred(i)
        mean = MVK[,,1]
        scale = MVK[,,2]
        krigvar = MVK[,,3]
      } else{
        for(j in 1:n.XX){
          # This is surprisingly slow, maybe there is a better way to do prediction
          # lagp_fit = lapply(1:n.pc,function(i) laGP::newGP(
          #   X = X[[i]][nn.indx[[i]][j,],,drop=F],
          #   Z = Z[i,nn.indx[[i]][j,]],
          #   d = 1,
          #   g = g))
          # lagp_pred = lapply(1:n.pc,function(i) laGP::predGP(lagp_fit[[i]],XX[[i]][j,,drop=F],lite=T))
          # Adding a few far away points can be helpful
          # Xadd = lapply(1:n.pc, function(i) sample(1:nrow(X[[i]]),10,F))
          Xadd = NULL

          for(i in 1:n.pc){
            lagp_fit = GP_fit_isotropic(rbind(X[[i]][nn.indx[[i]][j,],,drop=F],X[[i]][Xadd[[i]],]),
                                             d=1,g[i],
                                             c(Z[i,nn.indx[[i]][j,]],Z[i,Xadd[[i]]]),
                                             lite=F)
            lagp_pred = GP_predict(lagp_fit,
                                        rbind(X[[i]][nn.indx[[i]][j,],,drop=F],X[[i]][Xadd[[i]],]),
                                        XX[[i]][j,,drop=F],
                                        c(Z[i,nn.indx[[i]][j,]],Z[i,Xadd[[i]]]),
                                        predvar=T)
            mean[i,j] = lagp_pred$mean
            scale[i,j] = lagp_pred$s2
            krigvar[i,j] = lagp_pred$krigvar
          }
        }
      }
    }
  }
  # small nuggets can result in negative variances from laGP.
  if(predvar)
    scale[scale<0] = 0

  if(sample & predvar){
    sample = (rt(prod(dim(mean)),end) * sqrt(scale)) + mean
  } else{
    sample = mean
  }
  if(predvar){
    return(list(mean=mean,scale=scale,var=scale*end/(end-2),krigvar=krigvar,sample=sample,df=end,nn.indx=nn.indx))
  } else{
    return(list(mean=mean,sample=sample,df=end,nn.indx=nn.indx))
  }
}

garg = function(g, y){
  ## coerce inputs
  if(is.null(g)) g <- list()
  else if(is.numeric(g)) g <- list(start=g)
  if(!is.list(g)) stop("g should be a list or NULL")

  ## check mle
  if(is.null(g$mle)) g$mle <- FALSE
  if(length(g$mle) != 1 || !is.logical(g$mle))
    stop("g$mle should be a scalar logical")

  ## check for starting and max values
  if(is.null(g$start) || (g$mle && (is.null(g$max) || is.null(g$ab) || is.na(g$ab[2]))))
    r2s <- (y - mean(y))^2

  ## check for starting value
  # added that g$start >= sqrt(.Machine$double.eps) because that's the min, and there will be an error if start<min
  if(is.null(g$start)) g$start <- min(g$max-1e-12,max(sqrt(.Machine$double.eps),as.numeric(quantile(r2s, p=0.025))))
  ## check for max value
  if(is.null(g$max)) {
    if(g$mle) g$max <- max(r2s)
    else g$max <- max(g$start)
  }

  ## check for dmin
  if(is.null(g$min)) g$min <- sqrt(.Machine$double.eps)

  ## check for priors
  if(!g$mle) g$ab <- c(0,0)
  else {
    if(is.null(g$ab)) g$ab <- c(3/2, NA)
    if(is.na(g$ab[2])) {
      s2max <- mean(r2s)
      g$ab[2] <- laGP:::Igamma.inv(g$ab[1], 0.95*gamma(g$ab[1]), 1, 0)/s2max
    }
  }

  ## now check the validity of the values, and return
  # laGP:::check.arg(g)
  return(g)
}
