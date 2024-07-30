
# Returns predictions at X.pred.orig using the calibration parameters theta. To do this
# the emulator must be 'fit' at theta.
predict_w = function(flagp,X.pred.orig=NULL,theta=NULL,end=50,sample=F,n.samples=1,predvar=T)
{
  n.pc = flagp$basis$sim$n.pc
  y=NULL
  if(!is.null(X.pred.orig)){
    n.x.pred = nrow(X.pred.orig)
  } else{
    n.x.pred = 1
  }

  if(!is.null(theta)){
    X = transform_xt(X.sim = flagp$XT.data$sim$X$orig,
                     T.sim = flagp$XT.data$sim$T$orig,
                     X.obs = X.pred.orig,
                     T.obs = matrix(rep(theta,n.x.pred),nrow=n.x.pred,byrow=T))

    X = get_SC_inputs(flagp$lengthscales,X,n.pc)

    XX = lapply(1:n.pc,function(i) cbind(X$X.obs[[i]],X$T.obs[[i]]))
    w = aGPsep_SC_mv(X=X$XT.sim,
                     Z=flagp$basis$sim$V.t,
                     XX=XX,
                     start=start,
                     end=end,g=flagp$lengthscales$g,predvar=predvar)
  } else{
    X = transform_xt(X.sim = flagp$XT.data$sim$X$orig,
                     X.obs = X.pred.orig)

    X = get_SC_inputs(flagp$lengthscales,X,n.pc)
    XX = lapply(1:n.pc,function(i) X$X.obs[[i]])
    X = X$XT.sim

    ### DEV
    # Idea: screen variables by lengthscale: if we can remove unimportant inputs, nn calcs will be faster
    # ls.cuttoff = 10
    # for(i in 1:n.pc){
    #   keep.inputs = which(flagp$lengthscales$X[[i]]<ls.cuttoff)
    #   if(length(keep.inputs)==0){
    #     # none of the lengthscales were greater than the cuttoff, try just keeping the best input
    #     keep.inputs = which.max(flagp$lengthscales$X[[i]])
    #   }
    #   XX[[i]] = XX[[i]][,keep.inputs,drop=F]
    #   X[[i]] = X[[i]][,keep.inputs,drop=F]
    # }
    ### END DEV

    w = aGPsep_SC_mv(X=X,
                     Z=flagp$basis$sim$V.t,
                     XX=XX,
                     start=start,
                     end=end,bias=flagp$bias,
                     g=flagp$lengthscales$g,predvar=predvar)
  }
  if(sample){
    w$sample = t(sapply(1:n.samples, function (i) rt(prod(dim(w$mean)),w$df) * sqrt(w$var) + w$mean))
    dim(w$sample) = c(n.samples,dim(w$mean))
  } else{
    w$sample=w$mean
  }
  return(w)
}

mv_delta_predict = function(X.pred.orig,delta,flagp,sample=F,n.samples=1, start=6, end=50)
{
  n.pc = ncol(flagp$basis$obs$D)
  v = NULL
  v$df = delta$v$df
  # standardize X.pred.orig
  X.pred.std = unit_xform(X.pred.orig,X.min=flagp$XT.data$sim$X$min,X.range=flagp$XT.data$sim$X$range)$trans

  if(delta$method=='newGP'){
    # build full GPs
    delta.GPs = lapply(1:n.pc, function(k) laGP::newGPsep(X=flagp$XT.data$obs$X$trans,
                                                          Z=delta$V.t[k,],
                                                          d=delta$mle[[k]]$d, g=delta$mle[[k]]$g))
    # predict from GPs
    pred = lapply(1:n.pc, function(k) laGP::predGPsep(delta.GPs[[k]],X.pred.std,lite=T))
    for(i in 1:n.pc){
      v$mean = rbind(v$mean,pred[[i]]$mean)
      v$var = rbind(v$var,pred[[i]]$s2)
    }
  } else if(delta$method=='lagp'){
    # laGP model:
    # X - training inputs are observed X's
    # Z - training outputs are from the delta object which used residuals y - eta and basis D
    # XX - are the prediction locations scaled appropriately
    pred = aGPsep_SC_mv(X=flagp$XT.data$obs$X$trans,
                        Z=delta$V.t,
                        XX=X.pred.std,
                        start=start,
                        end=end,bias=T,
                        g=delta$mle$g)
    v$mean = pred$mean
    v$var = pred$var
  } else if(delta$method=='homgp'){
    pred = lapply(1:n.pc, function(k) predict(object = delta$model[[k]], x=X.pred.std))
    for(i in 1:n.pc){
      v$mean = rbind(v$mean, pred[[i]]$mean)
      v$var = rbind(v$var, pred[[i]]$sd2)
    }
  } else if(delta$method=='rgasp'){
    pred = lapply(1:n.pc, function(k) RobustGaSP::predict(delta$model[[k]],X.pred.std))
    for(i in 1:n.pc){
      v$mean = rbind(v$mean,pred[[i]]$mean)
      v$var = rbind(v$var,pred[[i]]$sd^2)
    }
  } else if(delta$method=='mlegp'){

  }

  if(sample){
    if(delta$method%in%c('lagp','newGP')){
      v$sample = t(sapply(1:n.samples, function (i) (rt(prod(dim(v$mean)),v$df) * sqrt(v$var)) + v$mean))
    } else{
      v$sample = t(sapply(1:n.samples, function (i) (rnorm(prod(dim(v$mean))) * sqrt(v$var)) + v$mean))
    }
    dim(v$sample) = c(n.samples,dim(v$mean)); v$sample = drop(v$sample)
  } else{
    v$sample=v$mean
  }
  return(v)
}

#' @title FlaGP Prediction
#'
#' @description Prediction with \code{mcmc} or \code{map} object
#' @param flagp an flagp data object.
#' @param model an \code{mcmc} or \code{map} object.
#' @param X.pred.orig new X for prediction on original scale
#' @param n.samples number of samples from the predictive distribution
#' @param samp.ids which mcmc samples should be used for prediction
#' @param return.samples return samples along with mean and confidence interval
#' @param support predict on support of simulations or observations
#' @param end.eta neighborhood size for lagp prediction
#' @param lagp.delta if a discrepancy model is fit, should discrepancy predictions be made with laGP? Default is a full GP.
#' @param start.delta initial lagp neighborhood size for discrepancy. Only relevant if lagp.delta = T
#' @param end.delta final lagp neighborhood size for discrepancy. Only relevant if lagp.delta = T
#' @param return.eta return emulator predictions
#' @param return.delta return discrepenacy predictions
#' @param y return predictions in scaled y space, rather than basis space
#' @param native return y predictions on native scale
#' @param conf.int return 95% confidence interval for predictions
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
predict.flagp = function(flagp,model=NULL,X.pred.orig=NULL,n.samples=0,samp.ids=NULL,return.samples=F,support='obs',
                         end.eta=50,lagp.delta=F,start.delta=6,end.delta=50,return.eta=F,return.delta=F,
                         y=T, native=T, conf.int=F,
                         resid.error = F, predvar = T, verbose=T)
{
  if(!is.null(model)){
    if(class(model)[1] == 'mcmc'){
      pred = mcmc_predict(flagp,model,X.pred.orig,samp.ids,n.samples,return.samples,support,end.eta,start.delta,end.delta,
                          return.eta,return.delta,native,conf.int,resid.error,predvar)
    } else if(class(model)[1] == 'map'){
      pred = map_predict(flagp,model,X.pred.orig,n.samples,return.samples,support,end.eta,start.delta,end.delta,
                         y,native,conf.int,resid.error,predvar)
    } else{
      stop('model must be of class mcmc or map')
    }
  } else{
    if(verbose)
      cat('No calibration model, emulation only prediction.')
    if(is.null(X.pred.orig))
      stop('must give X.pred.orig')
    pred = em_only_predict(flagp,X.pred.orig,n.samples,return.samples,support,end.eta,y,native,conf.int,predvar)
  }

  return(pred)
}

map_predict = function(flagp,map,X.pred.orig,n.samples,return.samples,support,
                       end.eta,start.delta,end.delta,y,native,conf.int,
                       resid.error,predvar,verbose)
{
  returns = list()
  if(support=='obs'){
    B = flagp$basis$obs$B
    D = flagp$basis$obs$D
    ym = flagp$Y.data$obs$mean
    ysd = flagp$Y.data$obs$sd
    n.y = flagp$Y.data$obs$n.y
  } else{
    B = flagp$basis$sim$B
    D = flagp$basis$sim$D
    ym = flagp$Y.data$sim$mean
    ysd = flagp$Y.data$sim$sd
    n.y = flagp$Y.data$sim$n.y
  }
  n.pc = ncol(B)
  n.pc.delta = ncol(D)

  if(!predvar & n.samples>0){
    stop('n.samples>0, set predvar=T')
  }

  if(resid.error){
    sigma = map$ssq.hat*diag(1,n.y)
  } else{
    sigma = .Machine$double.eps * diag(n.y)
  }

  # n = 1 if no X model
  n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
  # transform_theta
  theta = map$theta.hat * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min

  # emulator predictions
  start.time = proc.time()
  w = predict_w(flagp,X.pred.orig,theta,end=end.eta,predvar=predvar)
  returns$time = proc.time() - start.time

  if(!flagp$bias & y){
    # unbiased prediction
    if(predvar & n.samples>1)
      returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))

    for(i in 1:n.pred){
      if(predvar & n.samples>1){
        returns$y.samp[,,i] = t(t(mvnfast::rmvn(n.samples,B%*%w$mean[,i,drop=F],sigma + B%*%diag(w$var[,i]*w$df/(w$df-2),nrow=n.pc)%*%t(B))) *
                                  ysd + ym)
      } else{
        returns$y.mean = B%*%w$mean[,i,drop=F] * ysd + ym
      }
    }
    if(n.samples>1)
      returns$y.mean = apply(returns$y.samp,c(2,3),mean)
    if(conf.int)
      returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))
    if(!return.samples)
      returns$y.samp = NULL

  } else if(flagp$bias){
    # biased prediction
    if(n.samples>1){
      eps = sqrt(.Machine$double.eps)*diag(n.y)
      eta.samp = array(0,dim=c(n.y,n.samples,n.pred))
      delta.samp = array(0,dim=c(n.y,n.samples,n.pred))

      start.time = proc.time()
      v = mv_delta_predict(X.pred.orig,map$delta,flagp,F,start=start.delta,end=end.delta)
      returns$time = returns$time + proc.time() - start.time

      for(j in 1:n.pred){
        eta.samp[,,j] = t(mvnfast::rmvn(n.samples,B%*%w$mean[,j,drop=F],B%*%diag(w$var[,j]*w$df/(w$df-2),nrow=n.pc)%*%t(B) + eps)) * ysd + ym
        if(map$delta$method%in%c('lagp','newgp')){
          delta.samp[,,j] = t(mvnfast::rmvn(n.samples,D%*%v$mean[,j,drop=F],D%*%diag(v$var[,j]*v$df/(v$df-2),nrow=n.pc.delta)%*%t(D) + sigma )) * ysd # variance is scale of t * df/(df-2)
        } else{
          delta.samp[,,j] = t(mvnfast::rmvn(n.samples,D%*%v$mean[,j,drop=F],D%*%diag(v$var[,j],nrow=n.pc.delta)%*%t(D) + sigma )) * ysd
        }
      }
      if(y){
        returns$y.samp = eta.samp + delta.samp
        returns$y.mean = apply(returns$y.samp,c(1,3),mean) # do i need to compute this or just add eta.mean to delta.mean?
        if(conf.int)
          returns$y.conf.int = apply(returns$y.samp,c(1,3),quantile,c(.025,.975))
      }

      returns$eta.mean = apply(eta.samp,c(1,3),mean)
      returns$delta.mean = apply(delta.samp,c(1,3),mean)
      if(!return.samples)
        returns$y.samp = NULL
    } else{
      eta = array(0,dim=c(n.y,n.pred))
      delta = array(0,dim=c(n.y,n.pred))

      start.time = proc.time()
      v = mv_delta_predict(X.pred.orig,map$delta,flagp,F,start=start.delta,end=end.delta)
      returns$time = returns$time + proc.time() - start.time

      for(j in 1:n.pred){
        eta[,j] = t(B%*%w$mean[,j,drop=F]) * ysd + ym
        delta[,j] = t(D%*%v$mean[,j,drop=F]) * ysd
      }
      returns$y.mean = eta + delta
      returns$eta.mean = eta
      returns$delta.mean = delta
    }
  }
  return(returns)
}

mcmc_predict = function(flagp ,mcmc, X.pred.orig, samp.ids, n.samples, return.samples, support,
                        end.eta, start.delta, end.delta, return.eta, return.delta,
                        native, conf.int,
                        resid.error, predvar)
{
  returns = list()
  if(support=='obs'){
    B = flagp$basis$obs$B
    D = flagp$basis$obs$D
    ym = flagp$Y.data$obs$mean
    ysd = flagp$Y.data$obs$sd
    n.y = flagp$Y.data$obs$n.y
  } else{
    B = flagp$basis$sim$B
    D = flagp$basis$sim$D
    ym = flagp$Y.data$sim$mean
    ysd = flagp$Y.data$sim$sd
    n.y = flagp$Y.data$sim$n.y
  }
  if(is.null(samp.ids)){
    samp.ids = seq(1,(mcmc$n.samples - mcmc$n.burn), length.out = n.samples)
  } else{
    n.samples = length(samp.ids)
  }
  n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
  t.pred = t(t(mcmc$t.samp[samp.ids,,drop=F]) * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min)
  ssq.samp = mcmc$ssq.samp[samp.ids]
  returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
  returns$eta.samp = array(0,dim=c(n.samples,n.y,n.pred))
  if(flagp$bias)
    returns$delta.samp = array(0,dim=c(n.samples,n.y,n.pred))
  start.time = proc.time()
  for(i in 1:n.samples){
    w = predict_w(flagp,X.pred.orig,t.pred[i,],sample=T,end=end.eta)
    returns$eta.samp[i,,] = B%*%drop(w$sample)

    if(flagp$bias){
      # Biased prediction add delta model
      v = FlaGP:::mv_delta_predict(X.pred.orig,mcmc$delta[[i]],flagp,sample=T,n.samples=1,start=start.eta,end=end.eta)
      returns$delta.samp[i,,] = D%*%drop(v$sample) # if we scaled y.resid before fitting v, we probably need to scale again here
      for(j in 1:n.pred){
        if(resid.error){
          # sigma noise is on standardized scale, so add noise before scaling back to native with ysd and ym
          returns$y.samp[i,,j] = (returns$eta.samp[i,,j,drop=F]+returns$delta.samp[i,,j,drop=F] + (rnorm(n.y) * sqrt(ssq.samp[i])))
        } else{
          returns$y.samp[i,,j] = (returns$eta.samp[i,,j,drop=F]+returns$delta.samp[i,,j,drop=F])
        }
        if(native)
          returns$y.samp[i,,j] = returns$y.samp[i,,j] * ysd + ym
      }
    } else{
      # Unbiased prediction eta only
      for(j in 1:n.pred){
        if(resid.error){
          returns$y.samp[i,,j] = (returns$eta.samp[i,,j,drop=F] + (rnorm(n.y) * sqrt(ssq.samp[i])))
        } else{
          returns$y.samp[i,,j] = returns$eta.samp[i,,j,drop=F]
        }
        if(native)
          returns$y.samp[i,,j] = returns$y.samp[i,,j] * ysd + ym
      }
    }
  }
  returns$pred.time = proc.time() - start.time
  returns$y.mean = apply(returns$y.samp,2:3,mean)
  if(conf.int)
    returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))

  if(return.eta){
    # convert to standard scale first
    if(native){
      for(i in 1:n.samples){
        returns$eta.samp[i,,] = returns$eta.samp[i,,] * ysd + ym
      }
    }
    returns$eta.mean = apply(returns$eta.samp,2:3,mean)
    returns$eta.conf.int = apply(returns$eta.samp,2:3,quantile,c(.025,.975))
  }
  if(flagp$bias & return.delta){
    # convert to standard scale first
    if(native){
      for(i in 1:n.samples){
        returns$delta.samp[i,,] = returns$delta.samp[i,,] * ysd
      }
    }

    returns$delta.mean = apply(returns$delta.samp,2:3,mean)
    returns$delta.conf.int = apply(returns$delta.samp,2:3,quantile,c(.025,.975))
  }
  if(!return.samples)
    returns$y.samp = NULL
    returns$eta.samp = NULL
    if(flagp$bias)
      returns$delta.samp = NULL

  return(returns)
}

em_only_predict = function(flagp, X.pred.orig, n.samples, return.samples, support, end.eta, y, native, conf.int, predvar)
{
  returns = list()
  if(is.null(dim(X.pred.orig))){
    X.pred.orig = matrix(X.pred.orig,ncol=1)
  }
  n.pred = nrow(X.pred.orig)
  # get predictive samples of w at X.pred.orig
  start.time = proc.time()
  w = predict_w(flagp,X.pred.orig,end = end.eta,sample = ifelse(n.samples>0,T,F),n.samples = n.samples, predvar=predvar)
  returns$pred.time = proc.time() - start.time
  returns$w = w
  # convert samples of w to samples of y on native scale
  if(y){
    if(conf.int){
      returns$y.samp = array(0,dim=c(n.samples,flagp$Y.data$sim$n.y,n.pred))
      for(i in 1:n.samples){
        returns$y.samp[i,,] = flagp$basis$sim$B %*% drop(w$sample[i,,])
        if(native){
          for(j in 1:n.pred){
            returns$y.samp[i,,j] = returns$y.samp[i,,j] * flagp$Y.data$sim$sd + flagp$Y.data$sim$mean
          }
        }
      }
      returns$y.mean = apply(returns$y.samp,2:3,mean)
      returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(.025,.975))
    } else{
      returns$y.mean = flagp$basis$sim$B%*%returns$w$mean
      if(native)
        returns$y.mean = returns$y.mean * flagp$Y.data$sim$sd + flagp$Y.data$sim$mean
    }
    if(!return.samples)
      returns$y.samp = NULL
  }
  return(returns)
}

