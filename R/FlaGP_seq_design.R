# selects the X over the entire input space which has the largest prediction variance
seq_design_max_var_y = function(model,Xcand01,n.pc=model$basis$sim$n.pc,end=50){

  # TO DO
  # If the set of candidate points is not changing, we can skip recomputing the candidate point
  # variances that we have already computed in the previous iteration for all points where the
  # added point did not effect the neighborhood which can speed up the process.

  # predict from the model at all candidate locations
  pred = predict(model,X.pred.orig = Xcand01, verbose = F, y.var=T, X01=T,n.pc=min(n.pc,model$basis$sim$n.pc),end.eta = end)
  predvar = colMeans(pred$y.var)
  # return the X which has the largest prediction variance
  id = which.max(predvar)
  X = Xcand01[id,,drop=F]
  # IDEA: make an acquisition function that has a space filling property so that we don't just
  # pick points near the boundary
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,X_new_orig=X_orig,predvar=predvar,
              max_var = predvar[which.max(predvar)], cand_id = id))
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
seq_design_maximin = function(model,Xcand01){

  # compute distances for candidate points to each model point
  D = plgp::distance(Xcand01,model$XT.data$sim$X$trans)
  # find distance to closest model X for each cand X
  mins = apply(D,1,min)
  # pick the cand X with the maximum closest distance
  id = which.max(mins)

  X = Xcand01[id,,drop=F]
  X_orig = t((t(X) * model$XT.data$sim$X$range) + model$XT.data$sim$X$min)
  return(list(X_new=X,X_new_orig=X_orig,cand_id = id))
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
    return(seq_design_max_var_y(model,Xcand01,n.pc=n.pc,end=end))
  } else if(method=='rand'){
    return(seq_design_rand(model,Xcand01))
  } else if(method=='maximin'){
    return(seq_design_maximin(model,Xcand01))
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
        pred = predict(model,X.pred.orig = Xint[recompute,,drop=F], verbose = F,end.eta=end,y.var = T,n.pc=n.pc,X01=T,resid.error = F, parallel = parallel)
        imse_current_x[recompute] = colMeans(pred$y.var)
      }
    } else{
      pred = predict(model,X.pred.orig = Xint, verbose = F,end.eta=end,y.var = T,n.pc=n.pc,X01=T,resid.error = F, parallel = parallel)
      imse_current_x = colMeans(pred$y.var)
    }
    imse_current = mean(imse_current_x)

    # for each candidate X, compute the reduction in integrated variance when adding the candidate X to the training set
    imse_new = numeric(n_cand)

    # IMSE depends on y, so we need to add the new y (which we don't know), so lets predict it
    # unfortunately this prediction may be poor, not sure if there's any way around this. So, adding a point based on
    # a predicted y, then replacing that predicted y with the true y is likely to lead to instability in IMSE.
    pred_cand = predict(model,X.pred.orig = Xcand01, verbose = F,end.eta=end, y.var=F,n.pc=n.pc,X01=T)

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
        pred_int = predict(model_tmp,X.pred.orig = Xint[int_points_to_do[[i]],,drop=F], verbose = F,end.eta=end,y.var=T,n.pc=n.pc,X01=T, parallel = parallel)
        imse_cand_x[[i]][int_points_to_do[[i]]] = colMeans(pred_int$y.var)
      }
      imse_new[i] = mean(imse_cand_x[[i]])
    }

    d_imse = imse_current - imse_new
    which_cand = which.max(d_imse[d_imse>0])
    if(length(which_cand)==0){
      warning('No point reduced IMSE')
      return()
    }

    list(X_new = Xcand[which_cand,],
         cand_id = which_cand,
         imse_curr = imse_current,
         imse_new = imse_current - d_imse[which_cand],
         XintPredVars = imse_cand_x[[which_cand]],
         effected = int_points_to_do[[which_cand]])
  }
}
