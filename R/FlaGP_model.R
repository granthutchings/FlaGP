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
    delta = FlaGP:::fit_delta(eta$y.resid,flagp$basis$obs$D,flagp$XT.data,sample=sample,delta.method=delta.method,start=start.delta,end=end.delta,ssq.prior.params = ssq.prior.params,map=map)
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
# wrapper function for fit model that has theta and ssq contained in the first argument
fit_model_map = function(param,flagp,
                         end.eta=50,
                         delta.method='newGP',start.delta=6,end.delta=50,
                         theta.prior='beta',theta.prior.params=c(2,2),
                         ssq.prior='hcauchy',ssq.prior.params=c(.5)){
  theta = param[1:flagp$XT.data$p.t]
  ssq = param[flagp$XT.data$p.t+1]
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

  return(list(w=w,y.resid=y.resid))#,ssq.hat=ssq.hat))
}
# In mcmc we has previously recalculated the likelihood at the current theta to improve mixing.
# This is cumbersome when all we really needed was to resample the emulator and discrepancy model
resample_w = function(w,y.obs.trans,B.obs,ssq.prior.params){
  # resample w
  # w$sample = (rnorm(prod(dim(w$mean))) * sqrt(w$var)) + w$mean
  w$sample = rt(prod(dim(w$mean)),w$df) * sqrt(w$var) + w$mean
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
    #v$var = list()
    for(k in 1:n.pc){
      if(nrow(XT.data$obs$X$trans)>=5){
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
      # I suspect that these ranges and priors are bad for small n
      # cmle <- laGP::mleGPsep(GPs[[k]],
      #                        param='d',
      #                        tmin = d$min,
      #                        tmax = d$max,
      #                        ab=d$ab)
      cmle <- laGP::jmleGPsep(GPs[[k]],
                        drange=c(d$min, 10*d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      mle[[k]] = cmle
      pred = laGP::predGPsep(GPs[[k]], XX = XT.data$obs$X$trans, lite=T)

      # v$mean = rbind(v$mean, pred$mean * V.t.sd + V.t.mean)
      # v$var = rbind(v$var, pred$s2 / V.t.sd^2)

      #v$var[[k]] = pred$Sigma
      pred$s2[pred$s2<0]=0
      v$mean = rbind(v$mean, pred$mean)
      v$var = rbind(v$var, pred$s2)
      laGP::deleteGPsep(GPs[[k]])
    }
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
    # mlegp
    returns$model = vector(mode='list',length=n.pc)
    for(k in 1:n.pc){
      invisible(capture.output(returns$model[[k]] <- mlegp::mlegp(XT.data$obs$X$trans, V.t[k,],simplex.ntries = 1,verbose = 0, min.nugget = sqrt(.Machine$double.eps))))
      pred = predict(returns$model[[k]], se.fit = T)

      v$mean = rbind(v$mean, t(pred$fit))
      v$var = rbind(v$var, t(pred$se.fit)^2)
    }
  } else if(delta.method=='homgp'){
    # homgp
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
      v$sample = (rt(prod(dim(v$mean)),v$df) * sqrt(v$var)) + v$mean
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
resample_v = function(v,y.resid,delta.method){##ssq.prior.params,delta.method){
  # resample v
  if(delta.method %in% c('lagp','newGP')){
    v$sample = (rt(prod(dim(v$mean)),v$df) * sqrt(v$var)) + v$mean
  } else{
    v$sample = (rnorm(prod(dim(v$mean))) * sqrt(v$var)) + v$mean
  }
  # recompute residuals
  y.resid = y.resid - w_to_y(v$sample,D)# + v$y.resid.mean)

  # conjugate posterior
  # rtr = apply(y.resid,2,function(x) t(x) %*% x)
  # n.s2 = prod(dim(y.resid))
  # ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)

  return(list(v=v,y.resid=y.resid))#,ssq.hat=ssq.hat))
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
# var: matrix (n.pc x n) of prediction variances
## Calls:
# aGPsep from laGP package
# w_to_y to transform w -> y
## Called by:
# mv_calib.bias: biased calibration
# mv_calib.nobias: unbiased calibration
# mv.em.pred: laGP emulator prediction function at new inputs
}
aGPsep_SC_mv = function(X, Z, XX, g, start=6, end=50, bias=F, sample=F, predvar=T){

  n.pc = nrow(Z)
  n.XX = ifelse(bias,nrow(XX),nrow(XX[[1]]))
  # p.t = ifelse(bias,ncol(XX),ncol(XX[[1]]))
  mean = array(dim=c(n.pc,n.XX))
  var = array(dim=c(n.pc,n.XX))

  if(bias){
    # if we use laGP for the bias, we cannot use stretched an compressed inputs anymore so default to alc criterion
    lagp_fit = lapply(1:n.pc,function(i) laGP::aGPsep(X = X,
                                                Z = Z[i,],
                                                XX = XX,
                                                start = start,
                                                end = min(end,ncol(Z)-1),
                                                method = 'alc',
                                                verb=0))
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
                                d=1,g,lite=T)
          tKchol = t(GP$Kchol)
          # cross-cov: K[end+1,1:end]
          # train-cov: K[1:end,1:end]
          # traindata: Z[i,nn.indx[[i]][j,]]
          # Kt %*% KiY
          mean[i,j] = tKchol[end+1,1:end]%*%forwardsolve(tKchol[1:end,1:end],Z[i,nn.indx[[i]][j,]])
        }
      }
    } else{
      for(j in 1:n.XX){
        # This is surprisingly slow, maybe there is a better way to do prediction
        # lagp_fit = lapply(1:n.pc,function(i) laGP::newGP(
        #   X = X[[i]][nn.indx[[i]][j,],,drop=F],
        #   Z = Z[i,nn.indx[[i]][j,]],
        #   d = 1,
        #   g = g))
        # lagp_pred = lapply(1:n.pc,function(i) laGP::predGP(lagp_fit[[i]],XX[[i]][j,,drop=F],lite=T))
        for(i in 1:n.pc){
          lagp_fit = GP_fit_isotropic(X[[i]][nn.indx[[i]][j,],,drop=F],
                                      d=1,g,
                                      Z[i,nn.indx[[i]][j,]],
                                      lite=F)
          lagp_pred = GP_predict(lagp_fit,
                                 X[[i]][nn.indx[[i]][j,],,drop=F],
                                 XX[[i]][j,,drop=F],
                                 Z[i,nn.indx[[i]][j,]],
                                 predvar=T)
          mean[i,j] = lagp_pred$mean
          var[i,j] = lagp_pred$s2
        }
      }
    }
  }
  # small nuggets can result in negative variances from laGP.
  if(predvar)
    var[var<0] = 0

  if(sample & predvar){
    # sampling via cholesky where chol = sqrt(var)
    # sample = (rnorm(prod(dim(mean))) * sqrt(var)) + mean
    sample = (rt(prod(dim(mean)),end) * sqrt(var)) + mean
  } else{
    sample = mean
  }
  if(predvar){
    return(list(mean=mean,var=var,sample=sample,df=end))
  } else{
    return(list(mean=mean,sample=sample,df=end))
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
  if(is.null(g$start)) g$start <- max(sqrt(.Machine$double.eps),as.numeric(quantile(r2s, p=0.025)))
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
  laGP:::check.arg(g)
  return(g)
}
