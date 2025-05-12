
# w,v are needed if sample=F
compute_ll = function(theta,ssq,eta,delta,sample,flagp,
                      theta.prior,theta.prior.params,
                      ssq.prior,ssq.prior.params){
  ll = list()

  if(!is.null(delta)){
    # ssq.hat = delta$ssq.hat
    y.resid = delta$y.resid
  } else{
    # ssq.hat = eta$ssq.hat
    y.resid = eta$y.resid
  }
  if(sample){
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of y.resid
    ldetS = flagp$Y.data$obs$n.y*log(2*pi*ssq)
    # Don't compute S inverse directly - bad for huge d_y
    ySinvty = apply(y.resid^2,2,sum)/ssq
    ll = -.5*(ldetS + ySinvty)
  } else{
    # can't do this b/c cholesky of d_y x d_y
    # ll = numeric(ncol(y.resid))
    # mu = numeric(nrow(y.resid))
    # for(i in 1:n){
    #   ll[i] = mvnfast::dmvn(y.resid[,i],mu,flagp$basis$obs$B%*%diag(w.scale.adjust*eta$w$var[,i])%*%t(flagp$basis$obs$B) +
    #                           flagp$basis$obs$D%*%diag(v.scale.adjust*delta$v$var[,i])%*%t(flagp$basis$obs$D) +
    #                           diag(ssq,length(mu)),log = T)
    # }

    lambda = 1/ssq
    B.hat = flagp$precomp$BDtBDinvtBD%*%y.resid
    # yLLHmaty = diag(t(y.resid)%*%flagp$precomp$LLHmat%*%y.resid) # this is a bottleneck for large d_y
    yLLHmaty = colSums((flagp$precomp$LLHmat %*% y.resid) * y.resid) # this is the fast way of getting the diagonal, but still a slowdown for big d_y

    # need to make a variance adjustment for student-t
    w.scale.adjust = eta$w$df / (eta$w$df-2)
    if(flagp$bias){
      if(delta$method%in%c('lagp','newGP')){
        v.scale.adjust = delta$v$df / (delta$v$df-2)
        if(v.scale.adjust<=1) # small number of observations
          v.scale.adjust = 1
      } else{
        v.scale.adjust = 1
      }
    } else{
      v.scale.adjust = NULL
    }
    SigB.hat = lapply(1:flagp$Y.data$n, function(i) as.matrix(diag(c(w.scale.adjust*eta$w$var[,i],v.scale.adjust*delta$v$var[,i]),nrow=length(c(eta$w$var[,i],delta$v$var[,i]))) + ssq*flagp$precomp$BDtBDinv))
    # l_beta_hat = sapply(1:flagp$Y.data$n, function(i) mvtnorm::dmvnorm(x=as.numeric(B.hat[,i]),sigma = SigB.hat, log = T))
    mu = rep(0,nrow(B.hat))
    l_beta_hat = sapply(1:flagp$Y.data$n, function(i) mvnfast::dmvn(X=as.numeric(B.hat[,i]),mu=mu,sigma = SigB.hat[[i]], log = T))
    ll = sapply(1:flagp$Y.data$n, function(i) -.5*((flagp$Y.data$obs$n.y-flagp$precomp$rankBD)*log(2*pi*ssq) + flagp$precomp$ldetBDtBD +
                                                     lambda*yLLHmaty[i]) + l_beta_hat[i])
  }

  if(theta.prior=='beta'){
    p.theta = sum(dbeta(theta,theta.prior.params[1],theta.prior.params[2],log=T))
  } else if(theta.prior=='unif'){
    p.theta = sum(dunif(theta,0,1,log=T))
  }
  if(ssq.prior=='invgamma'){
    p.ssq = invgamma::dinvgamma(ssq,ssq.prior.params[1],ssq.prior.params[2],log=T)
  } else{
    # half cauchy
    p.ssq = LaplacesDemon::dhalfcauchy(ssq,ssq.prior.params[1],log=T)
  }
  return(sum(unlist(ll)) + p.theta + p.ssq)
}
