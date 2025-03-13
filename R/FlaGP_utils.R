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
