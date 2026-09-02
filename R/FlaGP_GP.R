GP_fit_isotropic = function(X,d,g,y=NULL,lite=F,XX=NULL)
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
  if(!lite & !is.null(XX)){
    returns$predXX = GP_predict(returns,X,XX,y,predvar=T)
  }
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


# drops a point in X, replacing it with Xnew, and returns the updated predictive variance at XX
GP_var_NN_student_t_update <- function(
    X, y, d, g, Xnew, XX, GP = NULL, dist = "L2", replace = TRUE
) {
  m <- nrow(X)
  if (m < 1) stop("X must have at least one row.")

  ## Helper: distance from a set to a single point
  dist_to_point <- function(M, vec, dist) {
    if (dist == "L2") laGP::distance(M, vec) else as.numeric(proxy::dist(M, vec, "manhattan"))
  }

  ## CHOLESKY solve helper
  solve_sym_pos <- function(M, v) {
    L <- chol(M)
    backsolve(L, forwardsolve(t(L), v))
  }

  ## ---------- REPLACE-ONE (drop + add) ----------
  if (isTRUE(replace)) {
    if (m <= 3) {  # df would be <= 0; fall back to add-only
      replace <- FALSE
    } else {
      # choose index to drop: farthest from XX
      Dxx <- dist_to_point(X, XX, dist)
      id_drop <- which.max(Dxx)

      R <- X[-id_drop, , drop = FALSE]                  # remaining (m-1) neighbors
      B <- GP_get_K(R, R, d, 1, g)                      # square block with nugget on diag
      L <- chol(B)
      solve_B <- function(v) backsolve(L, forwardsolve(t(L), v))

      # cross-covariances (no nugget on cross terms)
      u <- t(GP_get_K(R, XX,   d, 1, 0))                   # k(R, x)
      w <- solve_B(u)

      bprime <- t(GP_get_K(R, Xnew, d, 1, 0))            # k(R, x')
      qprime <- solve_B(bprime)
      sprime <- as.numeric(1 + g - crossprod(bprime, qprime))  # c' - b'^T B^{-1} b'
      tprime <- as.numeric(GP_get_K(Xnew, XX, d, 1, 0))        # k(x', x)
      gprime <- as.numeric(crossprod(bprime, w) - tprime)

      # unit-scale kriging variance after replace-one (predicting y, so k(x,x)=1+g)
      s2_base <- as.numeric(1 + g - (crossprod(u, w) + (gprime^2) / sprime))

      # Student-t expected scale after replacement is based on R (size m-1): phi_B/(m-3)
      yR <- y[-id_drop]
      zR <- solve_B(yR)
      phi_B <- as.numeric(crossprod(yR, zR))
      df <- (m - 1) - 2  # = m - 3 for zero-mean
      if (df <= 0) stop("Non-positive df in replace-one branch; set replace = FALSE.")
      s2_t <- (phi_B / df) * s2_base
      return(s2_t)
    }
  }

  ## ---------- ADD-ONLY (no drop) ----------
  # We reach here either because replace=FALSE or due to df.
  # Build A on the full current neighbor set D = X
  A <- GP_get_K(X, X, d, 1, g)
  L <- chol(A)
  solve_A <- function(v) backsolve(L, forwardsolve(t(L), v))

  # current unit-scale kriging variance at x: s2_old = (1+g) - u^T A^{-1} u
  u <- t(GP_get_K(X, XX, d, 1, 0))                    # k(D, x)
  Au <- solve_A(u)
  s2_old <- as.numeric((1 + g) - crossprod(u, Au))

  # insertion terms for x'
  b  <- t(GP_get_K(X, Xnew, d, 1, 0))                 # k(D, x')
  Ab <- solve_A(b)
  s  <- as.numeric(1 + g - crossprod(b, Ab))       # var_old(y(x'))
  t  <- as.numeric(GP_get_K(Xnew, XX, d, 1, 0))    # k(x', x)
  g_ins <- as.numeric(t - crossprod(b, Au))        # cov_old(y(x), y(x'))

  # unit-scale kriging variance after adding x'
  s2_base_new <- s2_old - (g_ins^2) / s

  # expected Student-t scale after adding equals c(D) = phi_A/(m-2)
  if (m <= 2) stop("Need at least 3 neighbors (m > 2) for Student-t variance in add-only.")
  Ay <- solve_A(y)
  phi_A <- as.numeric(crossprod(y, Ay))
  df_add <- m - 2                                   # zero-mean; with p-dim trend use m - p - 2
  s2_t_new <- (phi_A / df_add) * s2_base_new

  return(s2_t_new)
}
