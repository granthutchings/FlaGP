GP_fit_isotropic = function(X,d,g,y=NULL,lite=F)
{
  n = nrow(X)
  dx = ncol(X)

  out = .C("covar_symm_R",
           col = as.integer(dx),
           X = as.double(t(X)),
           n = as.integer(n),
           d = as.double(d),
           g = as.double(g),
           K = as.double(matrix(0,n,n)),
           PACKAGE = "FlaGP")
  dim(out$K) = c(n,n)
  returns = list(K = out$K,
                 Kchol = chol(out$K))
  if(!lite){
    if(is.null(y))
      stop('Must give y if lite=F')
    returns$Ki = chol2inv(returns$Kchol)
    returns$phi = t(y)%*%returns$Ki%*%y
  }
  returns$d = d
  returns$g = g
  returns$n = n
  returns$dx = dx
  return(returns)
}
# Given the cholesky of the covariance matrix draw a random sample
# from the mean zero with stderr 'sd'
GP_sample = function(gp_fit,sd)
{
  sd^2*drop(rnorm(nrow(gp_fit$Kchol)) %*% gp_fit$Kchol)
}

# With covariance K defined on points xy predict from the GP
# at XX using training points y
# XX must be a single point only
GP_predict = function(gp_fit,X,XX,y=NULL,predvar=F)
{
  # cross-covariance of training locations and prediction locations
  n1 = nrow(XX)
  n2 = nrow(X)
  dx = ncol(X)
  out = .C("covar_R",
           col = as.integer(dx),
           X1 = as.double(t(XX)),
           n1 = as.integer(n1),
           X2 = as.double(t(X)),
           n2 = as.integer(n2),
           d = as.double(gp_fit$d),
           K = as.double(matrix(0,n1,n2)),
           PACKAGE = "FlaGP")
  returns=list(K = out$K)
  dim(returns$K) = c(n1,n2)
  if(!is.null(y)){
    # mean prediction
    ktKi = returns$K%*%gp_fit$Ki
    returns$mean=ktKi%*%y

    if(predvar){
      phidf = gp_fit$phi/n2
      ktKik = tcrossprod(ktKi,returns$K)
      returns$krigvar = (1 + gp_fit$g)-ktKik # K_xx_xx = 1 + g
      returns$s2 = max(gp_fit$g,phidf*returns$krigvar)

    }
  }
  if(is.null(y) & predvar){
    stop('Can only give prediction variance if y is given.')
    # returns$s2 = diag(matrix(1+gp_fit$g) - returns$K %*% tcrossprod(gp_fit$Ki,returns$K))
  }
  return(returns)
}

GP_get_K = function(X,XX=NULL,d,s2,g){
  if(!is.null(XX)){
    n1 = nrow(XX)
    n2 = nrow(X)
    dx = ncol(X)
    out = .C("covar_R",
             col = as.integer(dx),
             X1 = as.double(t(XX)),
             n1 = as.integer(n1),
             X2 = as.double(t(X)),
             n2 = as.integer(n2),
             d = as.double(d),
             K = as.double(matrix(0,n1,n2)),
             PACKAGE = "FlaGP")
    dim(out$K) = c(n1,n2)
    out$K = out$K * s2
  } else{
    n = nrow(X)
    dx = ncol(X)
    out = .C("covar_symm_R",
             col = as.integer(dx),
             X = as.double(t(X)),
             n = as.integer(n),
             d = as.double(d),
             g = 0,
             K = as.double(matrix(0,n,n)),
             PACKAGE = "FlaGP")
    dim(out$K) = c(n,n)
    out$K = out$K * s2
    out$K = out$K + diag(g,n)
  }
  return(out$K)
}
