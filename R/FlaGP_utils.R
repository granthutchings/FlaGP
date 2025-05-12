makePositiveDefinite = function (m, tol)
{
  if (!is.matrix(m))
    m = as.matrix(m)
  d = dim(m)[1]
  if (dim(m)[2] != d)
    stop("Input matrix is not square!")
  es = eigen(m)
  esv = es$values
  if (missing(tol))
    tol = d * max(abs(esv)) * .Machine$double.eps
  delta = 2 * tol
  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
  return(m + dm)
}

invlogit = function (x)
{
  InvLogit <- 1/{
    1 + exp(-x)
  }
  return(InvLogit)
}

# this is the jacobian
dinvlogit = function(x,log=F){
  y = invlogit(x)
  z = y * (1 - y)
  if(log){
    return(log(z))
  } else{
    return(z)
  }
}

logit = function (p)
{
  if ({
    any(p < 0)
  } || {
    any(p > 1)
  })
    stop("p must be in [0,1].")
  Logit <- log(p/{
    1 - p
  })
  return(Logit)
}
