
# Returns predictions at X.pred.orig using the calibration parameters theta. To do this
# the emulator must be 'fit' at theta.
predict_w = function(flagp,X.pred.orig=NULL,theta=NULL,end=50,sample=F,n.samples=1,w.var=T,n.pc,parallel=F)
{
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
                     end=end,g=flagp$lengthscales$g,predvar=w.var)
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
                     Z=flagp$basis$sim$V.t[1:n.pc,,drop=F],
                     XX=XX,
                     start=start,
                     end=end,bias=flagp$bias,
                     g=flagp$lengthscales$g[1:n.pc],predvar=w.var,parallel=parallel)
  }
  if(sample){
    w$sample = t(sapply(1:n.samples, function (i) rt(prod(dim(w$mean)),w$df) * sqrt(w$scale) + w$mean))
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
      v$scale = rbind(v$scale,pred[[i]]$s2*(v$df-2)/v$df)
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
    v$scale = pred$scale
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
      v$sample = t(sapply(1:n.samples, function (i) (rt(prod(dim(v$mean)),v$df) * sqrt(v$scale)) + v$mean))
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
#' @param y.conf.int return 95% confidence interval for predictions
#' @param resid.error add residual standard error to predictions
#' @param w.var compute prediction variance in PC space
#' @param y.var compute prediction variance in y space
#' @param verbose print some information
#' @param alpha confidence level for confidence intervals
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
predict.flagp = function(flagp,model=NULL,X.pred.orig=NULL,n.samples=100,samp.ids=NULL,return.samples=F,support='obs',
                         end.eta=50,lagp.delta=F,start.delta=6,end.delta=50,return.eta=F,return.delta=F,
                         y=T, native=T, y.conf.int= T, joint = F,
                         resid.error = F, w.var = T, y.var = T, y.samp = T, verbose=T, alpha=.05, X01=F,
                         n.pc=flagp$basis$sim$n.pc, parallel = F, make.cluster = F)
{
  if(parallel & make.cluster){
    cl = parallel::makeCluster(min(flagp$basis$sim$n.pc,parallel::detectCores()-1))
    doParallel::registerDoParallel(cl)
  }
  if(!is.null(model)){
    if(class(model)[1] == 'mcmc'){
      if(verbose)
        cat('MCMC object given. Drawing predictive samples from posterior distribution.')
      # pred = mcmc_predict(flagp,model,X.pred.orig,samp.ids,n.samples,return.samples,support,end.eta,start.delta,end.delta,
      #                     return.eta,return.delta,native,y.conf.int,resid.error,w.var,y.var)
      pred = mcmc_predict_joint(flagp, model, X.pred.orig, n.samples, samp.ids,
                                support, native, joint,
                                end.eta, start.delta, end.delta,
                                return.eta, return.delta,
                                resid.error, y.var, y.samp, y.conf.int, alpha)
    } else if(class(model)[1] == 'map'){
      if(is.null(model$Cov)){
        if(verbose)
          cat('MAP object without covariance given. Drawing samples from predictive distribution at MLE.')
        pred = map_predict(flagp,model,X.pred.orig,n.samples,return.samples,support,end.eta,start.delta,end.delta,
                           y,native,y.conf.int,resid.error,w.var,y.var,alpha,return.eta,return.delta,verbose)
      } else{
        # we have the covariance matrix from the optimization, so take samples and do mcmc predict
        if(verbose)
          cat('MAP object with covariance given. Drawing samples from covariance matrix and making predictions.')
        pred = sample_predict(flagp, model, X.pred.orig, n.samples, return.samples, support,
                               end.eta, start.delta, end.delta, return.eta, return.delta,
                               native, y.conf.int,
                               resid.error, w.var,y.var)
      }

    } else{
      stop('model must be of class mcmc or map')
    }
  } else{
    if(verbose)
      cat('No calibration model, emulation only prediction.')
    if(is.null(X.pred.orig))
      stop('must give X.pred.orig')
    if(X01)
      X.pred.orig = t((t(X.pred.orig) * flagp$XT.data$sim$X$range) + flagp$XT.data$sim$X$min)
    pred = em_only_predict(flagp,X.pred.orig,n.samples,support,end.eta,y,native,return.samples,y.conf.int,w.var,y.var,alpha,n.pc,parallel)
  }

  return(pred)
}

map_predict = function(flagp,map,X.pred.orig,n.samples,return.samples,support,
                       end.eta,start.delta,end.delta,y,native,y.conf.int,
                       resid.error,w.var,y.var,alpha,return.eta,return.delta,verbose)
{
  get_yvar = function(B,wvar,ysd,s2,native=T,bias=F,D=NULL,vvar=NULL){
    diag_Sigma_w = diag(B%*%tcrossprod(diag(wvar),B))
    if(bias){
      diag_Sigma_v = diag(D%*%tcrossprod(diag(vvar,nrow=ncol(D)),D))
    } else{
      diag_Sigma_v = diag_Sigma_w * 0
    }

    if(!native)
      ysd = 1
    ysd^2*(s2 + diag_Sigma_w + diag_Sigma_v)
  }
  get_ysamp = function(B,wmean,wvar,s2,ysd,ym,n.samples,native=T,bias=F,D=NULL,vvar=NULL){
    Sigma_w = B%*%tcrossprod(diag(wvar,nrow=n.pc),B)
    if(bias){
      Sigma_v = D%*%tcrossprod(diag(vvar,nrow=ncol(D)),D)
    } else{
      Sigma_v = Sigma_w * 0
    }
    if(!native){
      ysd = 1
      ym = 0
    }

    diag(Sigma_w) = diag(Sigma_w) + s2  # add error variance
    t(t(mvnfast::rmvn(n.samples,B%*%wmean,Sigma_w+Sigma_v)) * ysd + ym)
  }
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

  if(resid.error){
    s2 = map$ssq.hat
  } else{
    s2 = sqrt(.Machine$double.eps)
  }
  if(!native){
    ysd = 1
    ym = 1
  }

  # n = 1 if no X model
  # and make a dummy x
  if(is.null(X.pred.orig)){
    n.pred = 1
    if(flagp$bias)
      X.pred.orig = matrix(.5,ncol=1)
  }

  # transform_theta
  theta = map$theta.hat * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min

  # emulator predictions
  start.time = proc.time()[3]
  w = predict_w(flagp,X.pred.orig,theta,end=end.eta,w.var=w.var,n.pc=n.pc)
  returns$pred.time = proc.time()[3] - start.time

  if(!flagp$bias & y){
    # mean
    returns$y.mean = B%*%w$mean * ysd + ym

    # variance, confidence interval, samples
    if(y.var | y.conf.int){
      returns$y.var = array(0,dim=dim(returns$y.mean))
      if(n.samples>1){
        returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
      }
      for(j in 1:n.pred){
        returns$y.var[,j] = get_yvar(B,w$var[,j],ysd,s2)
        if(n.samples>1){
          returns$y.samp[,,j] = get_ysamp(B,w$mean[,j,drop=F],w$var[,j],
                                          s2,ysd,ym,n.samples)
        }
      }
    }
  } else if(flagp$bias){
    # biased prediction
    start.time = proc.time()
    v = mv_delta_predict(X.pred.orig,map$delta,flagp,F,start=start.delta,end=end.delta)
    returns$time = returns$time + proc.time() - start.time

    # mean
    eta = B%*%w$mean * ysd + ym
    delta = D%*%v$mean * ysd
    returns$y.mean = eta + delta

    # variance, confidence interval, samples
    if(y.var | y.conf.int){
      returns$y.var = array(0,dim=dim(returns$y.mean))
      if(n.samples>1){
        returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
      }
      for(j in 1:n.pred){
        returns$y.var[,j] = get_yvar(B,w$var[,j],ysd,s2,bias = T,D = D,vvar = v$var[,j])
        if(n.samples>1){
          returns$y.samp[,,j] = get_ysamp(B,w$mean[,j,drop=F],w$var[,j],
                                          s2,ysd,ym,n.samples,bias=T,D=D,vvar=v$var[,j])
        }
      }
    }
  }
  if(y.conf.int){
    n.y = flagp$Y.data$sim$n.y
    returns$y.conf.int = array(dim=c(2,n.y,n.pred))
    returns$y.conf.int[1,,] = qnorm(alpha/2,returns$y.mean,sqrt(returns$y.var))
    returns$y.conf.int[2,,] = qnorm(1-(alpha/2),returns$y.mean,sqrt(returns$y.var))
  }
  if(return.eta)
    returns$eta = eta
  if(return.delta)
    returns$delta = delta
  return(returns)
}

mcmc_predict_joint = function(flagp, mcmc, X.pred.orig=NULL, n.samples=length(mcmc$ssq.samp), samp.ids=as.integer(seq(1,length(mcmc$ssq.samp),n.samples)),
                              support='sim', native=T, joint = F,
                              end.eta=min(flagp$num$m-1,50), start.delta = NULL, end.delta = NULL,
                              return.eta = F, return.delta = F,
                              resid.error = F, y.var = T, y.samp = T, y.conf.int = T, alpha = .05){
  get_ymean = function(B,wmean,ysd,ym,bias=F,D=NULL,vmean=NULL){
    mean = B%*%wmean
    if(!is.null(D) & !is.null(vmean)){
      mean = mean + D%*%vmean
    }
    return(ysd * mean + ym)
  }
  get_yvar = function(B,wvar,bias=F,D=NULL,vvar=NULL){
    Sigma_w = B%*%tcrossprod(diag(wvar),B)
    if(bias){
      Sigma_v = D%*%tcrossprod(diag(vvar,nrow=ncol(D)),D)
    } else{
      Sigma_v = 0
    }
    return(list(Sigma_w = Sigma_w, Sigma_v = Sigma_v))
  }
  get_ysamp = function(B,wmean,wvar,s2,ysd,ym,n.samples,bias=F,D=NULL,vvar=NULL){
    Sigma_w = B%*%tcrossprod(diag(wvar,nrow=n.pc),B)
    if(bias){
      Sigma_v = D%*%tcrossprod(diag(vvar,nrow=ncol(D)),D)
    } else{
      Sigma_v = Sigma_w * 0
    }
    diag(Sigma_w) = diag(Sigma_w) + s2  # add error variance
    t(t(mvnfast::rmvn(n.samples,B%*%wmean,Sigma_w+Sigma_v)) * ysd + ym)
  }

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
  bias = flagp$bias

  if(!native){
    ym = 0
    ysd = 1
  }
  if(is.null(samp.ids)){
    if(n.samples == 0)
      n.samples = 1
    samp.ids = as.integer(seq(1,length(mcmc$ssq.samp), length.out = n.samples))
  } else{
    n.samples = length(samp.ids)
  }
  n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
  t.pred = t(t(mcmc$t.samp[samp.ids,,drop=F]) * flagp$XT.data$sim$T$range + flagp$XT.data$sim$T$min)
  ssq.samp = mcmc$ssq.samp[samp.ids]

  returns$y.mean.samp = array(0,dim=c(n.samples,n.y,n.pred))
  if(y.var)
    returns$y.var.samp = array(0,dim=c(n.samples,n.y,n.pred))
  if(y.samp)
    returns$y.samp = array(0,dim=c(n.samples,n.y,n.pred))

  ptm = proc.time()[3]
  for(i in 1:n.samples){
    # get w
    w = predict_w(flagp,X.pred.orig,t.pred[i,],sample=F,w.var=T,end=end.eta,n.pc=n.pc)
    # get v
    if(flagp$bias){
      v = FlaGP:::mv_delta_predict(X.pred.orig,mcmc$delta[[i]],flagp,sample=F,start=start.delta,end=end.delta)
    } else{
      v = NULL
    }
    if(resid.error){
      Sigma_y = diag(ssq.samp[i],n.y)
    } else{
      Sigma_y = 0
    }
    # get mean
    returns$y.mean.samp[i,,] = get_ymean(B,w$mean,ysd,ym,bias,D,v$mean)

    if(y.var | y.samp){
      for(j in 1:n.pred){
        # get variance
        var = get_yvar(B,w$var[,j],bias,D,v$var[,j])
        Sigma = ysd^2*(var$Sigma_w + var$Sigma_v + Sigma_y)
        if(joint){
          CholSigma = chol(Sigma)
        } else{
          CholSigma = diag(sqrt(diag(Sigma)))
        }
        returns$y.var.samp[i,,j] = diag(Sigma)
        # get a sample
        returns$y.samp[i,,j] = mvnfast::rmvn(1,returns$y.mean.samp[i,,j],sigma = CholSigma,isChol = T)
      }
    }
  }
  returns$pred.time = ptm - proc.time()[3]
  returns$y.mean = apply(returns$y.samp,2:3,mean)
  returns$y.var = apply(returns$y.samp,2:3,var)
  if(y.conf.int){
    returns$y.conf.int = apply(returns$y.samp,2:3,quantile,c(alpha/2,1-alpha/2))
  }
  return(returns)
}

mcmc_predict = function(flagp ,mcmc, X.pred.orig, samp.ids, n.samples, return.samples, support,
                        end.eta, start.delta, end.delta, return.eta, return.delta,
                        native, y.conf.int,
                        resid.error, w.var,y.var){
  flagp = list(flagp)
  returns = list()
  for(k in 1:length(flagp)){
    returns[[k]] = list()
    if(support=='obs'){
      B = flagp[[k]]$basis$obs$B
      D = flagp[[k]]$basis$obs$D
      ym = flagp[[k]]$Y.data$obs$mean
      ysd = flagp[[k]]$Y.data$obs$sd
      n.y = flagp[[k]]$Y.data$obs$n.y
    } else{
      B = flagp[[k]]$basis$sim$B
      D = flagp[[k]]$basis$sim$D
      ym = flagp[[k]]$Y.data$sim$mean
      ysd = flagp[[k]]$Y.data$sim$sd
      n.y = flagp[[k]]$Y.data$sim$n.y
    }
    if(!native){
      ym = 0
      ysd = 1
    }
    if(is.null(samp.ids)){
      if(n.samples == 0)
        n.samples = 1
      samp.ids = as.integer(seq(1,length(mcmc$ssq.samp), length.out = n.samples))
    } else{
      n.samples = length(samp.ids)
    }
    n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
    t.pred = t(t(mcmc$t.samp[samp.ids,,drop=F]) * flagp[[k]]$XT.data$sim$T$range + flagp[[k]]$XT.data$sim$T$min)
    ssq.samp = mcmc$ssq.samp[samp.ids]
    returns[[k]]$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
    if(resid.error)
      returns[[k]]$y.samp.resid.error = array(0,dim=c(n.samples,n.y,n.pred))

    returns[[k]]$eta.samp = array(0,dim=c(n.samples,n.y,n.pred))
    if(flagp[[k]]$bias)
      returns[[k]]$delta.samp = array(0,dim=c(n.samples,n.y,n.pred))
    start.time = proc.time()
    for(i in 1:n.samples){
      # just get w, v predictions for timing
      w = predict_w(flagp[[k]],X.pred.orig,t.pred[i,],sample=T,end=end.eta,n.pc=flagp[[k]]$basis$sim$n.pc)
      returns[[k]]$eta.samp[i,,] = B%*%drop(w$sample)

      if(flagp[[k]]$bias){
        # Biased prediction add delta model
        v = FlaGP:::mv_delta_predict(X.pred.orig,mcmc$delta[[i]],flagp[[k]],sample=T,n.samples=1,start=start.delta,end=end.delta)
        returns[[k]]$delta.samp[i,,] = D%*%drop(v$sample) # if we scaled y.resid before fitting v, we probably need to scale again here
      }
    }
    returns[[k]]$pred.time = proc.time() - start.time

    for(i in 1:n.samples){
      # w = predict_w(flagp[[k]],X.pred.orig,t.pred[i,],sample=T,end=end.eta,n.pc=flagp[[k]]$basis$sim$n.pc)
      # returns[[k]]$eta.samp[i,,] = B%*%drop(w$sample)

      if(flagp[[k]]$bias){
        # Biased prediction add delta model
        # v = FlaGP:::mv_delta_predict(X.pred.orig,mcmc$delta[[i]],flagp[[k]],sample=T,n.samples=1,start=start.delta,end=end.delta)
        # returns[[k]]$delta.samp[i,,] = D%*%drop(v$sample) # if we scaled y.resid before fitting v, we probably need to scale again here
        for(j in 1:n.pred){
          if(resid.error){
            # sigma noise is on standardized scale, so add noise before scaling back to native with ysd and ym
            returns[[k]]$y.samp.resid.error[i,,j] = ((returns[[k]]$eta.samp[i,,j,drop=F]+returns[[k]]$delta.samp[i,,j,drop=F]) + (rnorm(n.y) * sqrt(ssq.samp[i]))) * ysd + ym
            returns[[k]]$y.samp[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F]+returns[[k]]$delta.samp[i,,j,drop=F]) * ysd + ym
          } else{
            returns[[k]]$y.samp[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F]+returns[[k]]$delta.samp[i,,j,drop=F]) * ysd + ym
          }
        }
      } else{
        # Unbiased prediction eta only
        for(j in 1:n.pred){
          if(resid.error){
            returns[[k]]$y.samp.resid.error[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F] + (rnorm(n.y) * sqrt(ssq.samp[i]))) * ysd + ym
            returns[[k]]$y.samp[i,,j] = returns[[k]]$eta.samp[i,,j,drop=F] * ysd + ym
          } else{
            returns[[k]]$y.samp[i,,j] = returns[[k]]$eta.samp[i,,j,drop=F] * ysd + ym
          }
        }
      }
    }
    returns[[k]]$y.mean = apply(returns[[k]]$y.samp,2:3,mean)
    if(y.conf.int)
      if(resid.error){
        returns[[k]]$y.conf.int = apply(returns[[k]]$y.samp.resid.error,2:3,quantile,c(.025,.975))
      } else{
        returns[[k]]$y.conf.int = apply(returns[[k]]$y.samp,2:3,quantile,c(.025,.975))
      }
    if(y.var)
      if(resid.error){
        returns[[k]]$y.var = apply(returns[[k]]$y.samp.resid.error,2:3,var)
      } else{
        returns[[k]]$y.var = apply(returns[[k]]$y.samp,2:3,var)
      }

    if(return.eta){
      # convert to standard scale first
      for(i in 1:n.samples){
        returns[[k]]$eta.samp[i,,] = returns[[k]]$eta.samp[i,,] * ysd + ym
      }
      returns[[k]]$eta.mean = apply(returns[[k]]$eta.samp,2:3,mean)
      returns[[k]]$eta.y.conf.int = apply(returns[[k]]$eta.samp,2:3,quantile,c(.025,.975))
    }

    if(flagp[[k]]$bias & return.delta){
      # convert to standard scale first
      for(i in 1:n.samples){
        returns[[k]]$delta.samp[i,,] = returns[[k]]$delta.samp[i,,] * ysd
      }
      returns[[k]]$delta.mean = apply(returns[[k]]$delta.samp,2:3,mean)
      returns[[k]]$delta.y.conf.int = apply(returns[[k]]$delta.samp,2:3,quantile,c(.025,.975))
    }

    if(!return.samples){
      returns[[k]]$y.samp = NULL
      returns[[k]]$y.samp.resid.error = NULL
    }
    if(!return.eta)
      returns[[k]]$eta.samp = NULL
    if(!return.delta)
      returns[[k]]$delta.samp = NULL
  }

  return(returns[[1]])
}

# calibrated prediction without having access to eta and delta samples from mcmc/map object
# this may be most useful when using MAP estimation with a Hessian derived covariance matrix
sample_predict = function(flagp, model, X.pred.orig, n.samples, return.samples, support,
                          end.eta, start.delta, end.delta, return.eta, return.delta,
                          native, y.conf.int,
                          resid.error, w.var,y.var){
  if(class(model)[1]=='mcmc'){
    stop('Only implemented for MAP models')
  }

  # samples in scaled space
  t_hat_scaled = LaplacesDemon::logit(model$theta.hat)
  ssq_hat_scaled = exp(model$ssq.hat)
  samples = mvnfast::rmvn(n.samples,c(t_hat_scaled,ssq_hat_scaled),model$Cov)
  # convert samples to native space
  samples[,1:flagp$num$p.t] = LaplacesDemon::invlogit(samples[,1:flagp$num$p.t])
  samples[,flagp$num$p.t+1] = exp(samples[,flagp$num$p.t+1])

  flagp = list(flagp)
  returns = list()
  for(k in 1:length(flagp)){
    returns[[k]] = list()
    if(support=='obs'){
      B = flagp[[k]]$basis$obs$B
      D = flagp[[k]]$basis$obs$D
      ym = flagp[[k]]$Y.data$obs$mean
      ysd = flagp[[k]]$Y.data$obs$sd
      n.y = flagp[[k]]$Y.data$obs$n.y
    } else{
      B = flagp[[k]]$basis$sim$B
      D = flagp[[k]]$basis$sim$D
      ym = flagp[[k]]$Y.data$sim$mean
      ysd = flagp[[k]]$Y.data$sim$sd
      n.y = flagp[[k]]$Y.data$sim$n.y
    }

    n.pred = ifelse(!is.null(X.pred.orig),nrow(X.pred.orig),1)
    t.pred = t(t(samples[,1:flagp[[k]]$XT.data$p.t]) * flagp[[k]]$XT.data$sim$T$range + flagp[[k]]$XT.data$sim$T$min)
    ssq.samp = samples[,flagp[[k]]$XT.data$p.t+1]
    returns[[k]]$y.samp = array(0,dim=c(n.samples,n.y,n.pred))
    returns[[k]]$eta.samp = array(0,dim=c(n.samples,n.y,n.pred))
    if(flagp[[k]]$bias)
      returns[[k]]$delta.samp = array(0,dim=c(n.samples,n.y,n.pred))
    start.time = proc.time()
    for(i in 1:n.samples){
      w = predict_w(flagp[[k]],X.pred.orig,t.pred[i,],sample=T,end=end.eta,n.pc=flagp[[k]]$basis$sim$n.pc)
      returns[[k]]$eta.samp[i,,] = B%*%drop(w$sample)

      if(flagp[[k]]$bias){
        # Biased prediction add delta model
        v = FlaGP:::mv_delta_predict(X.pred.orig,model$delta,flagp[[k]],sample=T,n.samples=1,start=start.delta,end=end.delta)
        returns[[k]]$delta.samp[i,,] = D%*%drop(v$sample) # if we scaled y.resid before fitting v, we probably need to scale again here
        for(j in 1:n.pred){
          if(resid.error){
            # sigma noise is on standardized scale, so add noise before scaling back to native with ysd and ym
            returns[[k]]$y.samp[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F]+returns[[k]]$delta.samp[i,,j,drop=F] + (rnorm(n.y) * sqrt(ssq.samp[i])))
          } else{
            returns[[k]]$y.samp[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F]+returns[[k]]$delta.samp[i,,j,drop=F])
          }
          if(native)
            returns[[k]]$y.samp[i,,j] = returns[[k]]$y.samp[i,,j] * ysd + ym
        }
      } else{
        # Unbiased prediction eta only
        for(j in 1:n.pred){
          if(resid.error){
            returns[[k]]$y.samp[i,,j] = (returns[[k]]$eta.samp[i,,j,drop=F] + (rnorm(n.y) * sqrt(ssq.samp[i])))
          } else{
            returns[[k]]$y.samp[i,,j] = returns[[k]]$eta.samp[i,,j,drop=F]
          }
          if(native)
            returns[[k]]$y.samp[i,,j] = returns[[k]]$y.samp[i,,j] * ysd + ym
        }
      }
    }
    returns[[k]]$pred.time = proc.time() - start.time
    returns[[k]]$y.mean = apply(returns[[k]]$y.samp,2:3,mean)
    if(y.conf.int)
      returns[[k]]$y.conf.int = apply(returns[[k]]$y.samp,2:3,quantile,c(.025,.975))
    if(y.var)
      returns[[k]]$y.var = apply(returns[[k]]$y.samp,2:3,var)

    if(return.eta){
      # convert to standard scale first
      if(native){
        for(i in 1:n.samples){
          returns[[k]]$eta.samp[i,,] = returns[[k]]$eta.samp[i,,] * ysd + ym
        }
      }
      returns[[k]]$eta.mean = apply(returns[[k]]$eta.samp,2:3,mean)
      returns[[k]]$eta.y.conf.int = apply(returns[[k]]$eta.samp,2:3,quantile,c(.025,.975))
    }
    if(flagp[[k]]$bias & return.delta){
      # convert to standard scale first
      if(native){
        for(i in 1:n.samples){
          returns[[k]]$delta.samp[i,,] = returns[[k]]$delta.samp[i,,] * ysd
        }
      }

      returns[[k]]$delta.mean = apply(returns[[k]]$delta.samp,2:3,mean)
      returns[[k]]$delta.y.conf.int = apply(returns[[k]]$delta.samp,2:3,quantile,c(.025,.975))
    }
    if(!return.samples)
      returns[[k]]$y.samp = NULL
    returns[[k]]$eta.samp = NULL
    if(flagp[[k]]$bias)
      returns[[k]]$delta.samp = NULL
  }

  return(returns[[1]])
}

em_only_predict = function(flagp, X.pred.orig, n.samples, support, end.eta, y, native, y.samp, y.conf.int, w.var, y.var, alpha, n.pc, parallel)
{
  returns = list()
  if(is.null(dim(X.pred.orig))){
    X.pred.orig = matrix(X.pred.orig,ncol=1)
  }
  if(y.conf.int & !y.var){
    y.var=T
  }

  n.pred = nrow(X.pred.orig)
  # get predictive samples of w at X.pred.orig
  start.time = proc.time()
  w = predict_w(flagp,X.pred.orig,end = end.eta,sample = ifelse(n.samples>0,T,F),n.samples = n.samples, w.var=w.var, n.pc=n.pc, parallel = parallel)
  returns$pred.time = proc.time() - start.time
  returns$w = w

  sd = 1
  mean = 0
  if(native){
    sd = flagp$Y.data$sim$sd
    mean = flagp$Y.data$sim$mean
  }

  # convert w to y
  if(y){
    returns$y.mean = flagp$basis$sim$B[,1:n.pc,drop=F]%*%w$mean * sd + mean
    if(y.var){
      returns$y.var = array(0,dim=dim(returns$y.mean))
      for(j in 1:n.pred){
        returns$y.var[,j] = diag(flagp$Y.data$sim$sd^2*flagp$basis$sim$B[,1:n.pc,drop=F]%*%tcrossprod(diag(w$var[1:n.pc,j],n.pc),flagp$basis$sim$B[,1:n.pc,drop=F]))
      }
    }
    # if y.samp
    if(y.samp){
      returns$y.samp = array(0,dim=c(n.samples,dim(returns$y.mean)))
      for(i in 1:n.samples){
        returns$y.samp[i,,] = flagp$basis$sim$B[,1:n.pc,drop=F] %*% returns$w$sample[i,,] * sd + mean
      }
    }
  }
  if(y.conf.int){
    n.y = flagp$Y.data$sim$n.y
    returns$y.conf.int = array(dim=c(2,n.y,n.pred))
    returns$y.conf.int[1,,] = qnorm(alpha/2,returns$y.mean,sqrt(returns$y.var))
    returns$y.conf.int[2,,] = qnorm(1-(alpha/2),returns$y.mean,sqrt(returns$y.var))
  }
  return(returns)
}

