# fit model at t=theta
fit_model = function(theta,flagp,
                     lite=T,sample=T,
                     end.eta=50,
                     lagp.delta=F,start.delta=6,end.delta=50,
                     negll=F,theta.prior='beta',theta.prior.params=c(2,2),
                     ssq.prior.params=c(1,1e-3),map=F){

  # Predict from emulator at [X.obs,theta]
  eta = FlaGP:::fit_eta(theta,flagp,sample,end.eta,ssq.prior.params,map=map)

  if(flagp$bias){
    delta = FlaGP:::fit_delta(eta$y.resid,flagp$basis$obs$D,flagp$XT.data,sample=sample,lagp=lagp.delta,start=start.delta,end=end.delta,ssq.prior.params = ssq.prior.params,map=map)
  } else{
    delta = NULL
  }

  ll = FlaGP:::compute_ll(theta,eta,delta,sample,flagp,theta.prior,theta.prior.params,ssq.prior.params)

  if(lite){
    return(ifelse(negll,-ll,ll))
  } else{
    if(flagp$bias){
      ssq.hat = delta$ssq.hat
      eta$ssq.hat = NULL; delta$ssq.hat = NULL
    } else{
      ssq.hat = eta$ssq.hat
      eta$ssq.hat = NULL
    }

    return(list(ll = ifelse(negll,-ll,ll),
                theta = theta,
                ssq.hat = ssq.hat,
                eta = eta,
                delta = delta))
  }
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
                   sample=sample)

  # residuals
  y.pred = w_to_y(w$sample,flagp$basis$obs$B)
  y.resid = flagp$Y.data$obs$trans - y.pred

  # conjugate posterior
  rtr = apply(y.resid,2,function(x) t(x) %*% x)
  n.s2 = prod(dim(flagp.data$Y.data$obs$trans))
  if(!map){
    ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)
  } else{
    # mean of inv gamma for MAP
    ssq.hat = (ssq.prior.params[2]+sum(rtr)/2)/(ssq.prior.params[1]+n.s2/2-1)
  }

  return(list(w=w,y.resid=y.resid,ssq.hat=ssq.hat))
}

# fit delta model to residuals using basis vectors in D
fit_delta = function(y.resid,D,XT.data,sample=T,lagp=F,start=6,end=50,ssq.prior.params,map){
  returns = list()
  v = NULL
  n.pc = ncol(D)
  V.t = get_basis(y.resid,B=D)$V.t
  if(lagp){ # use laGP for a faster bias model
    v = append(v,aGPsep_SC_mv(X=XT.data$obs$X$trans,
                              Z=V.t,
                              XX=XT.data$obs$X$trans,
                              start=start,
                              end=end,
                              bias=T))#,SC=F))
  } else{       # use full GP for the bias model
    GPs = vector(mode='list',length=n.pc)
    mle = vector(mode='list',length=n.pc)
    #v$var = list()
    for(k in 1:n.pc){
      d = laGP::darg(d=list(mle=T),X=XT.data$obs$X$trans)
      g = garg(g=list(mle=T),y=V.t[k,])
      # Our responses have changed w.r.t to inputs, so we can't not use d=1 anymore
      GPs[[k]] <- laGP::newGPsep(X=XT.data$obs$X$trans,
                             Z=V.t[k,],
                             d=d$start, g=g$start, dK=TRUE)
      cmle <- laGP::jmleGPsep(GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      mle[[k]] = cmle
      pred = laGP::predGPsep(GPs[[k]], XX = XT.data$obs$X$trans, lite=T)
      v$mean = rbind(v$mean, pred$mean)
      v$var = rbind(v$var, pred$s2)
      #v$var[[k]] = pred$Sigma
      laGP::deleteGPsep(GPs[[k]])
    }
    returns$mle=mle
  }
  # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
  if(sample){
    # sampling via cholesky where chol = sqrt(var)
    v$sample = (rnorm(prod(dim(v$mean))) * sqrt(v$var)) + v$mean
    #for(k in 1:n.pc){
    #  v$sample = rbind(v$sample,mvtnorm::rmvnorm(1,mean=v$mean[k,],sigma=v$var[[k]]))
    #}
  } else{
    v$sample = v$mean
  }
  returns$v = v
  returns$V.t = V.t
  returns$y.resid = y.resid - w_to_y(v$sample,D)

  rtr = apply(returns$y.resid,2,function(x) t(x) %*% x)
  n.s2 = prod(dim(flagp.data$Y.data$obs$trans))
  if(!map){
    returns$ssq.hat = invgamma::rinvgamma(1,ssq.prior.params[1]+n.s2/2,ssq.prior.params[2]+sum(rtr)/2)
  } else{
    # mean of inv gamma for MAP
    returns$ssq.hat = (ssq.prior.params[2]+sum(rtr)/2)/(ssq.prior.params[1]+n.s2/2-1)
  }
  returns$lagp = lagp
  return(returns)
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
aGPsep_SC_mv = function(X, Z, XX, start=6, end=50, g=1/10000, bias=F, sample=F, ann = T){

  n.pc = nrow(Z)
  n.XX = ifelse(bias,nrow(XX),nrow(XX[[1]]))
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
    if(!ann){
      lagp_fit = lapply(1:n.pc,function(i) laGP::aGP(
        X = X[[i]],
        Z = Z[i,],
        XX = XX[[i]],
        d = list(mle = FALSE, start = 1),
        g = g,
        end = min(end,ncol(Z)-1),
        method = 'nn',
        verb=0))

      for (i in 1:n.pc) {
        mean[i,] = lagp_fit[[i]]$mean
        var[i,] = lagp_fit[[i]]$var
      }
    } else{
      # compute neighbors using ANN -- for small n.pc (2) and large m (10000) lapply is faster than parallel, probably overhead with parallel
      nn.indx = lapply(1:n.pc,function(i) yaImpute::ann(X[[i]],XX[[i]],end,verbose=F)$knnIndexDist[,1:end,drop=F])
      # foreach can be faster depending on n.pc and the number of available cores. We have not tried to explicitly dermine where this tradeoff lives
      # nn.indx = foreach::foreach(i=1:n.pc) %dopar% yaImpute::ann(X[[i]],XX[[i]],end,verbose=F)$knnIndexDist[,1:end]

      for(j in 1:n.XX){
        lagp_fit = lapply(1:n.pc,function(i) laGP::newGP(
          X = X[[i]][nn.indx[[i]][j,],,drop=F],
          Z = Z[i,nn.indx[[i]][j,]],
          d = 1,
          g = g))
        lagp_pred = lapply(1:n.pc,function(i) laGP::predGP(lagp_fit[[i]],XX[[i]][j,,drop=F],lite=T))
        for(i in 1:n.pc){
          mean[i,j] = lagp_pred[[i]]$mean
          var[i,j] = lagp_pred[[i]]$s2
          laGP::deleteGP(lagp_fit[[i]])
        }
        # using newGP/predGP limits us to the exponential covariance. If we write these GP functions by hand we can use the matern
        # Ex.
        # nu = 7/2
        # C = fields::Matern(distances::distances(X[[i]][nn.indx[[i]][j,],]),smoothness = nu)
      }
    }
  }

  if(sample){
    # sampling via cholesky where chol = sqrt(var)
    sample = (rnorm(prod(dim(mean))) * sqrt(var)) + mean
  } else{
    sample = mean
  }
  return(list(lagp_fit=lagp_fit,mean=mean,var=var,sample=sample))
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
