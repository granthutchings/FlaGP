# selects the X over the entire input space which has the largest prediction variance
seq_design_max_var_y = function(model,Xcand01,n.pc=model$basis$sim$n.pc,end=50,parallel=T){

  # TO DO
  # If the set of candidate points is not changing, we can skip recomputing the candidate point
  # variances that we have already computed in the previous iteration for all points where the
  # added point did not effect the neighborhood which can speed up the process.

  # predict from the model at all candidate locations
  pred = predict(model,X.pred.orig = Xcand01, verbose = F, y.var=T, X01=T,
                 n.pc=min(n.pc,model$basis$sim$n.pc),NN = end, parallel=parallel)
  predvar = colMeans(pred$y.var)
  # return the X which has the largest prediction variance
  id = which.max(predvar)
  X = Xcand01[id,,drop=F]
  # IDEA: make an acquisition function that has a space filling property so that we don't just
  # pick points near the boundary
  # X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,predvar=predvar,
              max_var = predvar[id], cand_id = id))
}
# selects random X over the entire input space from a candidate set
seq_design_rand = function(model,Xcand01){
  n = nrow(Xcand01)
  id = sample(1:n,1)
  X = Xcand01[id,,drop=F]
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,X_new_orig=X_orig,cand_id = id))
}
# selects X over the entire input space from a candidate set to respect maximin with the previous set of points
seq_design_maximin = function(model,Xcand01,scaled=F,n.pc=model$basis$sim$n.pc,dist='L2'){

  if(scaled){
    if(n.pc>1){
      # do maximin w.r.t. lengthscale estimates, using a weighted distance based on eigen value
      XcandSCj = lapply(1:n.pc, function(j) FlaGP:::sc_inputs(Xcand01,model$lengthscales$XT[[j]]))
      if(dist=='L1'){
        Dj = lapply(1:n.pc, function(j) proxy::dist(XcandSCj[[j]],model$SC.inputs$XT.sim[[j]],'manhattan'))
      } else{
        Dj = lapply(1:n.pc, function(j) sqrt(plgp::distance(XcandSCj[[j]],model$SC.inputs$XT.sim[[j]])))
      }
      w = model$basis$sim$singular.values/sum(model$basis$sim$singular.values)
      D = matrix(0,nrow=nrow(Dj[[1]]),ncol=ncol(Dj[[1]]))
      for(j in 1:n.pc){
        D = D + w[j] * Dj[[j]]
      }
    } else{
      # do maximin w.r.t. 1st PC lengthscale estimates - you'd think this would be better
      XcandSC = FlaGP:::sc_inputs(Xcand01,model$lengthscales$XT[[1]])
      if(dist=='L1'){
        D = proxy::dist(XcandSC,model$SC.inputs$XT.sim[[1]],'manhattan')
      } else{
        D = sqrt(plgp::distance(XcandSC,model$SC.inputs$XT.sim[[1]]))
      }
    }
  } else{
    # compute distances for candidate points to each model point
    if(dist=="L1"){
      D = proxy::dist(Xcand01,model$XT.data$sim$X$trans,'manhattan')
    } else{
      D = sqrt(plgp::distance(Xcand01,model$XT.data$sim$X$trans))
    }
  }
  # find distance to closest model X for each cand X
  mins = apply(D,1,min)
  # pick the cand X with the maximum closest distance
  id = which.max(mins)

  X = Xcand01[id,,drop=F]
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,X_new_orig=X_orig,cand_id = id))
}
# ALC criterion using only the first GP
seq_design_alc_1 = function(model,Xcand,Xint,NN=50){

  nC = nrow(Xcand)
  nI = nrow(Xint)

  XcandS = sweep(Xcand,2,sqrt(model$lengthscales$XT[[1]]),'/')
  XintS = sweep(Xint,2,sqrt(model$lengthscales$XT[[1]]),'/')
  XtrainS = model$SC.inputs$XT.sim[[1]]
  g = model$lengthscales$g[1]
  jitter = 1e-8
  p = ncol(XtrainS)
  d = rep(1,p)

  # get nearest neighbors for candidate points
  if(NN<nrow(XtrainS)){
    nn_index = FNN::get.knnx(XtrainS,XcandS,NN)$nn.index
  } else{
    nn_index = matrix(rep(seq(1:nrow(XtrainS)),nC),nrow=nC,byrow = T)
  }
  m = ncol(nn_index)

  # for each candidate, compute the ALC over the integration set
  alc = numeric(nC)
  for (j in 1:nC) {
    idx <- nn_index[j, ]
    XNs  <- XtrainS[idx, , drop = FALSE]   # m x p (scaled)
    xps  <- XcandS[j, , drop = FALSE] # 1 x p (scaled)

    # --- Local correlation matrices (no sigma^2) ---
    # R(XN, XN): via squared distances in scaled space
    K = FlaGP:::GP_get_K(XNs,d = d,s2 = 1,g = g)

    # r(XN, x')
    Kxx = FlaGP:::GP_get_K(XNs,XX = xps,d = d,s2 = 1,g = g)

    # s_corr^2(x') = 1 - r_*^T K^{-1} r_*
    # s2_corr_xp <- max(1 - sum(r_star * t), 0.0)
    s2_corr_xp = max(1-(Kxx%*%chol2inv(chol(K))%*%t(Kxx)),0.0)
    denom <- s2_corr_xp + g

    # r(z, x') for all z in Zint
    # we should only be computing this for those x in Xint where adding x' would effect their neighborhood
    if(NN<(nrow(XtrainS)+1)){
      # if NN<#train+1(cand) only compute c(z,x') for those points in z that are effected
      nn_z = FNN::get.knnx(rbind(xps,XtrainS),XintS,NN)$nn.index
      int_points_to_do <- which(apply(nn_z == 1, 1, any)) # if 1 is in the nn set, include it
    } else{
      # if NN>=number of training data, compute cz for all
      int_points_to_do = 1:nrow(XintS)
    }
    if(length(int_points_to_do)>0){
      # 3) for each affected z, compute c(z, x') using the z-local system
      cz2_sum <- 0
      k_z = min(NN,nrow(XtrainS))
      XintSi = Xint[int_points_to_do,]
      # neighbors of z from training only (exclude x'); size = k_pred - 1
      idx_z <- FNN::get.knnx(XtrainS, XintSi, k = k_z)$nn.index
      for (ri in 1:length(int_points_to_do)) {
        z <- XintSi[ri, , drop = FALSE]

        XNz   <- XtrainS[idx_z[ri,], , drop = FALSE]

        # K(z): R(XNz, XNz) + gI
        K_z = FlaGP:::GP_get_K(XNz,XX = NULL,d = d,s2 = 1,g = g)
        L_z   <- chol(K_z)

        # vectors for z and x'
        r_z = t(FlaGP:::GP_get_K(XNz,XX = z,d = d,s2 = 1,g = g))
        r_zxp <- t(FlaGP:::GP_get_K(XNz,XX = xps,d = d,s2 = 1,g = g))

        # c(z, x') = r(z,x') - R(z,XNz) K(z)^{-1} r(XNz, x')
        # Here R(z, XNz) = r_z^T (length k_z)
        t1   <- forwardsolve(t(L_z), r_zxp)
        tvec <- backsolve(L_z, t1)
        cz   <- as.numeric( exp(-sum((z - xps)^2)) - crossprod(r_z, tvec) )

        cz2_sum <- cz2_sum + (cz * cz)
      }
      alc[j] <- cz2_sum / length(int_points_to_do) / denom
    }
  }

  id = which.max(alc)
  X = Xcand[id,,drop=F]
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,X_new_orig=X_orig,cand_id=id,alc=alc))
}

seq_design_alc_1_replace <- function(model, Xcand, Xint, NN = 50) {
  nC <- nrow(Xcand)
  nI <- nrow(Xint)

  ## ----- scaling to correlation space (like your original) -----
  XcandS  <- sweep(Xcand,  2, sqrt(model$lengthscales$XT[[1]]), "/")
  XintS   <- sweep(Xint,   2, sqrt(model$lengthscales$XT[[1]]), "/")
  XtrainS <- model$SC.inputs$XT.sim[[1]]
  g       <- model$lengthscales$g[1]
  p       <- ncol(XtrainS)
  d       <- rep(1, p)  # kernel uses pre-scaled inputs, so set d=1

  ## ----- precompute: for each candidate, a training NN set (not strictly needed for replacement) -----
  if (NN < nrow(XtrainS)) {
    nn_index_cand <- FNN::get.knnx(XtrainS, XcandS, k = NN)$nn.index
  } else {
    nn_index_cand <- matrix(rep(seq_len(nrow(XtrainS)), nC), nrow = nC, byrow = TRUE)
  }

  ## ALC scores
  alc <- numeric(nC)

  for (j in seq_len(nC)) {
    xps <- XcandS[j, , drop = FALSE]             # candidate (scaled)

    ## ----- find affected integration points: x' enters z's NN set (replacement) -----
    # we check neighbors among {x'} ∪ XtrainS, and select z where index 1 (x') is included
    if (NN <= nrow(XtrainS)) {
      nn_z <- FNN::get.knnx(rbind(xps, XtrainS), XintS, k = NN)$nn.index
      int_points_to_do <- which(apply(nn_z == 1, 1, any))
    } else {
      # if NN > n_train, x' will always enter (but we'll still DROP the farthest training neighbor)
      int_points_to_do <- seq_len(nI)
    }

    if (length(int_points_to_do) == 0) {
      alc[j] <- 0
      next
    }

    ## accumulate Δ_z over affected z
    delta_sum <- 0
    k_z <- min(NN, nrow(XtrainS))                # size of z's training NN set (before replacement)

    # Pre-fetch z's training NN indices (size k_z)
    idx_z_all <- FNN::get.knnx(XtrainS, XintS[int_points_to_do, , drop = FALSE], k = k_z)$nn.index

    for (ri in seq_along(int_points_to_do)) {
      z_idx <- int_points_to_do[ri]
      z     <- XintS[z_idx, , drop = FALSE]

      # z's neighbors from training only (size m = k_z)
      idx_z <- idx_z_all[ri, ]
      XNz   <- XtrainS[idx_z, , drop = FALSE]

      # identify the point to drop: the farthest neighbor from z within XNz
      dists <- rowSums((XNz - matrix(z, nrow = nrow(XNz), ncol = p, byrow = TRUE))^2)
      j_drop <- which.max(dists)
      dj     <- XNz[j_drop, , drop = FALSE]

      # R = XNz without the dropped point
      if (nrow(XNz) >= 2) {
        R  <- XNz[-j_drop, , drop = FALSE]
        # B = K(R,R) + gI
        B  <- FlaGP:::GP_get_K(R, d = d, s2 = 1, g = g)
        L  <- chol(B)
        # helper: solve(B, v)
        solve_B <- function(v) backsolve(L, forwardsolve(t(L), v))

        # u = k(R, z), common term; w = B^{-1} u
        u  <- t(FlaGP:::GP_get_K(R, XX = z,   d = d, s2 = 1, g = g))
        w  <- solve_B(u)

        # OLD term (with dropped neighbor d_j):
        b   <- t(FlaGP:::GP_get_K(R, XX = dj, d = d, s2 = 1, g = g))
        q   <- solve_B(b)
        s   <- as.numeric(1 + g - crossprod(b, q))               # c - b^T B^{-1} b
        # t  = k(d_j, z) in correlation space:
        t   <- as.numeric(exp(-sum((dj - z)^2)))
        g_sc <- as.numeric(crossprod(b, w) - t)                  # g = b^T w - t

        # NEW term (with candidate x'):
        bprime <- t(FlaGP:::GP_get_K(R, XX = xps, d = d, s2 = 1, g = g))
        qprime <- solve_B(bprime)
        sprime <- as.numeric(1 + g - crossprod(bprime, qprime))  # c' - b'^T B^{-1} b'
        tprime <- as.numeric(exp(-sum((xps - z)^2)))             # k(x', z)
        gprime <- as.numeric(crossprod(bprime, w) - tprime)      # g' = b'^T w - t'

        # variance reduction at z under replacement:
        delta_z <- (g_sc * g_sc) / s - (gprime * gprime) / sprime

      } else {
        # Degenerate m = 1 case: R is empty. Use closed form:
        t   <- as.numeric(exp(-sum((dj  - z)^2)))
        tpr <- as.numeric(exp(-sum((xps - z)^2)))
        s   <- 1 + g
        sprime <- 1 + g
        delta_z <- (t * t - tpr * tpr) / s
      }

      # numeric guard
      if (delta_z > 0) delta_sum <- delta_sum + delta_z
    }

    # ALC score for candidate j: average over affected z
    alc[j] <- delta_sum / length(int_points_to_do)
  }

  id <- which.max(alc)
  X_best <- Xcand[id, , drop = FALSE]
  X_best_orig <- t((t(X_best) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)

  list(
    X_new       = X_best,
    X_new_orig  = X_best_orig,
    cand_id     = id,
    alc         = alc
  )
}

seq_design_alc_t_rank_1_updates <- function(model, Xcand, Xint, NN, dist = "L2", s2_curr=rep(NA,nrow(Xint))) {

  Xcand = sweep(Xcand,2,sqrt(model$lengthscales$XT[[1]]),'/')
  Xint = sweep(Xint,2,sqrt(model$lengthscales$XT[[1]]),'/')
  X = model$SC.inputs$XT.sim[[1]]
  y = model$basis$sim$V.t[1,]
  g = model$lengthscales$g[1]

  stopifnot(is.matrix(X), is.matrix(Xcand), is.matrix(Xint))
  n <- nrow(X); p <- ncol(X)
  d = rep(1,p)

  stopifnot(ncol(Xcand) == p, ncol(Xint) == p)
  if (n < 1) stop("X must have at least one training point.")
  if (NN < 1) stop("NN must be >= 1.")

  # distance helper: distances from all rows of M to a single row vec
  dist_to_point <- function(M, vec=NULL, dist='L2') {
    if (dist == "L2") {
      laGP::distance(M, vec)
    } else if (dist == "L1") {
      as.numeric(proxy::dist(M, vec, "manhattan"))
    } else {
      stop("dist must be 'L2' or 'L1'.")
    }
  }

  # current Student-t variance at a single integration point using its k nearest neighbors
  current_t_var <- function(X, y, XX) {
    GP_fit_isotropic(X,d,g,y,XX=XX)
  }

  # --- Precompute neighbor sets, thresholds, and current variances over Xint ---
  k0 <- min(NN, n)                               # current neighborhood size actually used
  replace_mode <- (NN <= n)                      # TRUE: drop+add, FALSE: add-only
  n_int <- nrow(Xint)

  nn_idx_list <- vector("list", n_int)
  rmax <- numeric(n_int)                         # distance to the k0-th neighbor

  for (j in 1:length(s2_curr)) {
    dj <- dist_to_point(X, Xint[j, , drop = FALSE], dist)
    idx <- order(dj)[1:k0]
    nn_idx_list[[j]] <- idx
    rmax[j] <- if (replace_mode) dj[idx[k0]] else Inf  # if add-only, always include new point
    if(is.na(s2_curr[j])){
      Xnb <- X[idx, , drop = FALSE]
      ynb <- y[idx]
      s2_curr[j] <- current_t_var(Xnb, ynb, Xint[j, , drop = FALSE])$predXX$s2
    }
  }

  # --- Evaluate ALC gain for each candidate ---
  k_cand <- nrow(Xcand)
  imspe <- numeric(k_cand)

  for (i in seq_len(k_cand)) {
    # gain_i <- 0.0
    xnew <- Xcand[i, , drop = FALSE]
    s2_tmp = s2_curr
    for (j in seq_len(n_int)) {
      # Will xnew be used at Xint[j,] under our neighborhood policy?
      # - replace_mode: only if closer than current farthest neighbor
      # - add-only: always TRUE (rmax[j] = Inf)
      d_cand <- dist_to_point(xnew, Xint[j, , drop = FALSE], dist)  # scalar
      if (d_cand < rmax[j]) {
        idx <- nn_idx_list[[j]]
        Xnb <- X[idx, , drop = FALSE]
        ynb <- y[idx]

        s2_new <- GP_var_NN_student_t_update(
          X = Xnb, y = ynb, d = d, g = g,
          Xnew = xnew, XX = Xint[j, , drop = FALSE],
          GP = NULL, dist = dist, replace = replace_mode
        )
        s2_tmp[j] = s2_new
        # delta <- s2_curr[j] - as.numeric(s2_new)
        # if (delta > 0) gain_i <- gain_i + delta   # clamp negatives (numerics) to zero
      }
    }
    imspe[i] <- mean(s2_tmp)
  }

  best_idx <- which.min(imspe)
  X = Xcand[best_idx,,drop=F]
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  list(
    cand_id = best_idx,
    X_new = X,
    X_new_orig = X_orig,
    best_imspe = imspe[best_idx],
    imspes = imspe,
    s2_current = s2_curr,
    nn_indices = nn_idx_list,
    rmax = rmax,
    replace_mode = replace_mode,
    neighborhood_size_used = k0
  )
}


#' @title FlaGP sequential design
#'
#' @description selects the X which reduces IMSE over the entire domain when added to the training set
#' @param model FlaGP model to be augmented with a new simulation point
#' @param n_cand number of candidate points
#' @param n_int number of integration points
#' @param Xcand optional matrix of candidate points bounded between 0 and 1
#' @param Xint optional matrix of integration points bounded between 0 and 1
#' @param seed seed for reproducibility
#' @param end number of neighborhood points for prediction
#' @param n.pc number of PC's to use for design. Less can be much faster
#' @details Returns new X point
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
seq_design = function(model,n_cand=100,n_int=100,
                      Xcand = NULL,Xint = NULL,
                      XintPredVars = NULL,
                      seed=NULL, end=50, n.pc=-1, method='imse',
                      dist='L2',
                      parallel = T,
                      verbose=T){
  if(!is.null(seed))
    set.seed(seed)

  if(n.pc==-1 | n.pc>model$basis$sim$n.pc){
    # use whatever number of components the model has
    n.pc = model$basis$sim$n.pc
  }
  if(n.pc == 1)
    parallel = F
  # compute MSE over entire support
  XcandGiven=F
  if(!is.null(Xcand))
    XcandGiven = T
  XintGiven=F
  if(!is.null(Xint))
    XintGiven=T

  if(!XcandGiven){
    # if no Xcand given, use a design on [0,1]
    Xcand01 = lhs::maximinLHS(n_cand,model$num$p.x + model$num$p.t)
  } else{
    Xcand01 = Xcand
    n_cand = nrow(Xcand01)
  }
  if(method=='maxvar' | method == 'MaxVar'){
    return(seq_design_max_var_y(model,Xcand01,n.pc=n.pc,end=end,parallel=parallel))
  } else if(method=='rand'){
    return(seq_design_rand(model,Xcand01))
  } else if(method %in% c('maximin','maximin-scaled')){
    return(seq_design_maximin(model,Xcand01,scaled = ifelse(method=='maximin-scaled',T,F),n.pc=n.pc,dist=dist))
  } else if(method=='alc'){
    return(seq_design_alc_1(model,Xcand01,Xint,end))
  } else{
    if(!XintGiven){
      Xint = lhs::maximinLHS(n_int,model$num$p.x + model$num$p.t)
    } else{
      n_int = nrow(Xint)
    }

    # transform from 0,1 to native space for prediction
    if(model$num$p.t>0){
      Xcand = t((t(Xcand01) * model$XT.data$sim$XT$range) + model$XT.data$sim$XT$min)
    } else{
      Xcand = t((t(Xcand01) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
    }

    # scale integration, candidate, and training points points
    XintSC = lapply(1:model$basis$sim$n.pc, function(i) FlaGP:::sc_inputs(Xint,model$lengthscales$XT[[i]]))
    XcandSC = lapply(1:model$basis$sim$n.pc, function(i) FlaGP:::sc_inputs(Xcand01,model$lengthscales$XT[[i]]))
    XtrainSC = lapply(1:model$basis$sim$n.pc, function(i) FlaGP:::sc_inputs(cbind(model$XT.data$sim$X$trans,model$XT.data$sim$T$trans),model$lengthscales$XT[[i]]))

    end = min(end,model$num$m+1)

    # we can compute the variance of y in closed form from the variance of w
    # this means that we could also get the quantiles of y
    if(verbose)
      cat('\ncomputing current IMSE ...\n')

    if(!is.null(XintPredVars)){
      # some prediction variances have already been computed under this model, so don't recompute them
      imse_current_x = XintPredVars
      recompute = which(is.na(imse_current_x))
      if(length(recompute)>0){
        pred = predict(model,X.pred.orig = Xint[recompute,,drop=F], verbose = F,NN=end,y.var = T,n.pc=n.pc,X01=T,resid.error = F, parallel = parallel)
        imse_current_x[recompute] = colMeans(pred$y.var)
      }
    } else{
      pred = predict(model,X.pred.orig = Xint, verbose = F,NN=end,y.var = T,n.pc=n.pc,X01=T,resid.error = F, parallel = parallel)
      imse_current_x = colMeans(pred$y.var)
    }
    imse_current = mean(imse_current_x)

    # for each candidate X, compute the reduction in integrated variance when adding the candidate X to the training set
    imse_new = numeric(n_cand)

    # IMSE depends on y, so we need to add the new y (which we don't know), so lets predict it
    # unfortunately this prediction may be poor, not sure if there's any way around this. So, adding a point based on
    # a predicted y, then replacing that predicted y with the true y is likely to lead to instability in IMSE.
    pred_cand = predict(model,X.pred.orig = Xcand01, verbose = F,NN=end, y.var=F,n.pc=n.pc,X01=T, parallel = parallel)

    if(verbose){
      cat('searching candidate space ... \n')
      pb = txtProgressBar(min = 1, max = n_cand/2, initial = 1)
    }
    imse_cand_x = int_points_to_do = list()
    for(i in 1:n_cand){
      if(verbose)
        setTxtProgressBar(pb,i/2)

      # make copy of the model object that can be changed
      model_tmp = model


      # add candidate point and predicted y to model
      model_tmp = FlaGP::flagp_update(model_tmp,Xcand[i,,drop=F],NULL,pred_cand$y.mean[,i],refit=F)

      # determine which integration points contain the candidate point in their NN set over any of the PC's
      # we need to recompute all the prediction variances for these effected points
      int_points_to_do[[i]] <- tryCatch({
        if(n.pc > 1){
          nn.indx <- abind::abind(lapply(1:n.pc,
                                         function(j) {FNN::get.knnx(rbind(XcandSC[[j]][i,, drop = FALSE], XtrainSC[[j]]),XintSC[[j]],end)$nn.index}),along = 0)
          int_points_to_do[[i]] <- which(apply(nn.indx == 1, 2, any))
        } else{
          nn.indx <- FNN::get.knnx(rbind(XcandSC[[1]][i,, drop = FALSE], XtrainSC[[1]]),XintSC[[1]],end)$nn.index
          int_points_to_do[[i]] <- which(apply(nn.indx == 1, 1, any))
        }
        int_points_to_do[[i]]
      }, error = function(e){
        message("KNN error: ", e$message)
        return(0)
      })

      # compute new IMSE
      imse_cand_x[[i]] = imse_current_x
      if(length(int_points_to_do[[i]])>0){
        pred_int = predict(model_tmp,X.pred.orig = Xint[int_points_to_do[[i]],,drop=F], verbose = F,NN=end,y.var=T,n.pc=n.pc,X01=T, parallel = parallel)
        imse_cand_x[[i]][int_points_to_do[[i]]] = colMeans(pred_int$y.var)
      }
      imse_new[i] = mean(imse_cand_x[[i]])
    }

    d_imse = imse_current - imse_new
    which_cand = which.max(d_imse)
    if(length(which_cand)==0){
      warning(paste0('No point reduced IMSE, accepting point that increases imse by ',d_imse[which_cand],'\n'))
    }

    list(X_new = Xcand[which_cand,],
         cand_id = which_cand,
         imse_curr = imse_current,
         imse_new = imse_current - d_imse[which_cand],
         XintPredVars = imse_cand_x[[which_cand]],
         effected = int_points_to_do[[which_cand]])
  }
}
