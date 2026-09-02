library(ggplot2)
library(Rtools)
library(FlaGP)
setwd("~/FlaGP/examples/multivariate_sepia/sequential_design")
NN = 10
refitevery = 1
n.pc = 5
lagp = F
do_flagp = T
plot_during_seqd = T
save_results = F

##########################
method = 'maximin-scaled'
n.pc.sd = 1
##########################

if(method=='maximin-scaled')
  lagp = F
if(lagp){
  source('lagp_pred.R')
  start = 1
}


# notes
# seed 735184 created the figure in the paper right now
# seed 6112025 - names for the date, this one was run using 5 PC always rather than going from 3 to 5
# seed 6202025 - recreating 735184 without the change in n.pc with the addition of MAPE
seed = 6202025
save_file = paste0(method,'_results_seed_',seed,'_refit_',refitevery,'_PC_',n.pc.sd,'_smart_refit_flagp_only.RData')
if(lagp==T & do_flagp==F)
  save_file = paste0(method,'_results_seed_',seed,'_refit_',refitevery,'_PC_',n.pc.sd,'smart_refit.RData')

if(n.pc.sd == 'all')
  n.pc.sd = -1

parallel = T
if(parallel){
  cl = parallel::makeCluster(n.pc)
  doParallel::registerDoParallel(cl)
}

nugget = 1e-10
# load initial ensemble
D = reticulate::py_load_object('ensemble_sequential_design.pkl')
D$y_sim = D$y_sim + 10
par(mfrow=c(2,2))
matplot(t(D$y_sim),type='l')

N = nrow(D$t_sim)
set.seed(seed)
# initialize with something like an LHS
# initialize points using those close to edges and middle
init_points = as.matrix(expand.grid(c(0,1),c(0,1),c(0,1)))
init_points = rbind(init_points,c(.5,.5,.5))
m_init = nrow(init_points)
init_ids = c()
for(i in 1:nrow(init_points)){
  d = plgp::distance(init_points[i,,drop=F],D$t_sim)
  init_ids = c(init_ids,which.min(d))
}
cand_ids = sample((1:N)[-init_ids],500,F)
test_ids = (1:N)[-c(init_ids,cand_ids)]

D.init = list()
D.init$X = D$t_sim[init_ids,]
D.init$Y = t(D$y_sim[init_ids,])
D.init$ids = init_ids
matplot(D.init$Y,type='l')

D.test = list()
D.test$X = D$t_sim[test_ids,]
D.test$Y = t(D$y_sim[test_ids,])
D.test$ids = test_ids
matplot(D.test$Y,type='l')

D.cand = list()
D.cand$X = D$t_sim[cand_ids,]
D.cand$Y = t(D$y_sim[cand_ids,])
D.cand$ids = cand_ids
matplot(D.cand$Y,type='l')

X_int = lhs::maximinLHS(1000,3)
n_points_to_add = 30-m_init
imse = imse_flagp = rmse = rmse_flagp = mape = mape_flagp = numeric(n_points_to_add+1)

matplot(D.init$Y,type='l')
ls.subsample = 'all'
ls.m = 2
ls.K = 1
model = flagp(X.sim = D.init$X, Y.sim = D.init$Y, n.pc = n.pc, center=T, scale=T,
              ls.subsample = ls.subsample, ls.K = ls.K, ls.m = ls.m, transform_x = F,
              verbose=F, ls.nugget = nugget, nug.est = F, ls.parallel = F, make.cluster = F,seed = 11)
plot(model)
cumsum(model$basis$sim$pct.var)

if(lagp){
  lagp.model.init = list(yo=model$Y.data$sim$orig,yt=model$Y.data$sim$trans,
                         xo=cbind(model$XT.data$sim$X$orig,model$XT.data$sim$T$orig),
                         xt=cbind(model$XT.data$sim$X$trans,model$XT.data$sim$T$trans),
                         xmin = c(model$XT.data$sim$X$min,model$XT.data$sim$T$min), xrange = c(model$XT.data$sim$X$range,model$XT.data$sim$T$range),
                         B=model$basis$sim$B, w=model$basis$sim$V.t, ym=model$Y.data$sim$mean, ysd=model$Y.data$sim$sd,
                         num=list(p.x=model$num$p.x, p.t = model$num$p.t,m=model$num$m), n.pc = n.pc, do_transform=model$XT.data$transformed)
  pred_int = lagp.predict.PCA(lagp.model.init,X_int,start,min(lagp.model.init$num$m-1,NN),nug = nugget, nug_est = F)
  imse[1] = mean(pred_int$y.var)
  pred_test = lagp.predict.PCA(lagp.model.init,D.test$X,start,min(lagp.model.init$num$m-1,NN),nug = nugget, nug_est = F)
  mape[1] = Rtools:::get_mape(pred_test$y.mean,D.test$Y)
  rmse[1] = Rtools:::get_rmse(pred_test$y.mean,D.test$Y)
}

flagp.int.pred = predict(model,X.pred.orig = X_int, X01=T,y.var = T,NN = NN, resid.error = F)
imse_flagp[1] = mean(flagp.int.pred$y.var)
flagp.test.pred = predict(model,X.pred.orig = D.test$X, X01=T,y.var = T,NN = NN, resid.error = F)
rmse_flagp[1] = Rtools:::get_rmse(flagp.test.pred$y.mean,D.test$Y)
mape_flagp[1] = Rtools:::get_mape(flagp.test.pred$y.mean,D.test$Y)
lagp.time = flagp.time = flagp.design.time = numeric(n_points_to_add+1)

flagp.time[1] = model$time
par(mfrow=c(1,3))
set.seed(1)
# we dont have the data generating mechanism, so the candidate and integration sets must be part of the known ensemble
xtcand = list(D.cand$X)
xtint = X_int
xcandflagp = xtcand[[1]]
if(lagp){
  lagp.model = list(); lagp.model[[1]] = lagp.model.init
}
flagp.model = list();flagp.model[[1]] = model
pred.lagp = list()
force_refit = F

for(i in 1:n_points_to_add){
  cat(i,' ')
  m = m_init + i
  refit = ifelse(m%%refitevery==0,T,F)
  if(m==10)
    refit = F
  if(force_refit)
    refit = T

  ####################################################################################################
  # FlaGP
  if(do_flagp){
    ptm = proc.time()
    if(i==1)
      XintVars = NULL

    # find new point
    if(method=='alc'){
      xnewF = FlaGP:::seq_design_alc_t_rank_1_updates(model = flagp.model[[max(1,i-1)]],
                                                     Xcand = xcandflagp,
                                                     Xint = xtint,
                                                     NN = NN,dist = 'L2')
    } else{
      xnewF = FlaGP::seq_design(model = flagp.model[[max(1,i-1)]],
                                Xcand = xcandflagp,
                                Xint = xtint,
                                XintPredVars = XintVars,
                                end = NN,method=method,n.pc = n.pc.sd,
                                verbose=F, parallel = parallel)
    }


    flagp.design.time[i+1] = -(ptm - proc.time())[3]

    # remove the row from xcand and XcandVars
    XintVars = xnewF$XintPredVars
    xcandflagp = xcandflagp[-xnewF$cand_id,]

    # find y associated with new point
    idnew = Rtools:::row_match(matrix(xnewF$X_new,nrow=1),D.cand$X)
    yneworig = D.cand$Y[,idnew$rowid]

    # update model with new point
    flagp.model[[i]] = FlaGP::flagp_update(flagp.model[[max(1,i-1)]],matrix(xnewF$X_new,nrow=1),Ynew_orig = yneworig,refit = refit,
                                           ls.K=ls.K,ls.subsample=ls.subsample,ls.m=ls.m,seed=3*i,n.pc=n.pc,parallel=parallel,make.cluster=F)
    force_refit = F # always set force_refit to false because we either just tried refitting, or we don't want to refit yet
    # force_refit only becomes true if we choose to ignore a bad refit

    # flag candidate points that need variance recalculated
    if(refit | method == 'maxvar'){
      # the model was refit so we need to recalculate everything
      XintVars = NULL
    } else if(method == 'imse'){
      XintVars[xnewF$effected] = NA
    }
    flagp.time[i+1] = -(ptm - proc.time())[3]
  }

  ####################################################################################################
  ####################################################################################################
  if(lagp){
    ptm = proc.time()
    if(method=='IMSE' | method=='imse'){
      xnewLaGP = lagp_seq_design_imse(lagp.model[[max(1,i-1)]],Xcand = xtcand[[i]],Xint = xtint,
                                      start = start, end = min(lagp.model[[max(1,i-1)]]$num$m-1,NN),
                                      nugget = nugget, nug_est = F, verbose=F)
      idnew = Rtools:::row_match(matrix(xnewLaGP$X_new,nrow=1),D.cand$X)
      yneworig = D.cand$Y[,idnew$rowid]
      lagp.model[[i]] = lagp_update(lagp.model[[max(1,i-1)]], xneworig = xnewLaGP$X_new, yneworig = yneworig, refit = refit, n.pc = n.pc)
      lagp.time[i+1] = -(ptm - proc.time())[3]
      xtcand[[i+1]] = xtcand[[i]][-xnewLaGP$cand_id,]
    } else if(method=='maxvar' | method=='maximin'){
      if(method=='maximin'){
        id = FlaGP:::seq_design_maximin(lagp.model[[max(1,i-1)]],xtcand[[i]])$cand_id
      } else{
        pred.lagp[[i]] = lagp.predict.PCA(lagp.model[[max(1,i-1)]],xtcand[[i]],start,min(lagp.model[[max(1,i-1)]]$num$m-1,NN),
                                          nug = nugget,nug_est = F)
        id = which.max(colSums(pred.lagp[[i]]$y.var))
      }
      xtadd = xoadd = xtcand[[i]][id,,drop=F]
      idnew = Rtools:::row_match(matrix(xtadd,nrow=1),D.cand$X)
      yneworig = D.cand$Y[,idnew$rowid]
      lagp.model[[i]] = lagp_update(lagp.model[[max(1,i-1)]], xneworig = xoadd, yneworig = yneworig, refit = refit, n.pc = n.pc)
      lagp.time[i+1] = -(ptm - proc.time())[3]
      xtcand[[i+1]] = xtcand[[i]][-id,]
    }
  }

  ####################################################################################################

  # imse / rmse
  if(lagp){
    lagp.pred = lagp.predict.PCA(lagp.model[[i]],X_int,start = 6,end = min(lagp.model[[i]]$num$m-1,NN),nug = nugget, nug_est = F)
    imse[i+1] = mean(lagp.pred$y.var)
    lagp.pred = lagp.predict.PCA(lagp.model[[i]],D.test$X,start = 6,end = min(lagp.model[[i]]$num$m-1,NN),nug = nugget, nug_est = F)
    rmse[i+1] = Rtools:::get_rmse(lagp.pred$y.mean,D.test$Y)
    mape[i+1] = Rtools:::get_mape(lagp.pred$y.mean,D.test$Y)
  }
  if(do_flagp){
    flagp.pred = predict(flagp.model[[i]],X.pred.orig = X_int, X01=T,y.var = T,verbose = F,NN = NN, parallel = parallel)
    imse_flagp[i+1] = mean(flagp.pred$y.var)
    flagp.pred = predict(flagp.model[[i]],X.pred.orig = D.test$X, X01=T,y.var = T,verbose = F,NN = NN, parallel = parallel)
    rmse_flagp[i+1] = Rtools:::get_rmse(flagp.pred$y.mean,D.test$Y)
    mape_flagp[i+1] = Rtools:::get_mape(flagp.pred$y.mean,D.test$Y)
    cat('rmse new:',rmse_flagp[i+1],'\n')
    if(refit & rmse[i+1]>(1.1*rmse[i])){
      cat('bad refit, reverting to last fit, and trying refit again next iteration\n')
      # we refit, and it got worse, so don't refit
      flagp.model[[i]] = FlaGP::flagp_update(flagp.model[[max(1,i-1)]],matrix(xnewF$X_new,nrow=1),Ynew_orig = yneworig,refit = F,
                                             seed=k*i,n.pc=flagp.model[[max(1,i-1)]]$basis$sim$n.pc,parallel=parallel,make.cluster=F)
      flagp.pred = predict(flagp.model[[i]],X.pred.orig = X_int, X01=T,y.var = T,verbose = F,NN = NN, parallel = parallel)
      imse_flagp[i+1] = mean(flagp.pred$y.var)
      flagp.pred = predict(flagp.model[[i]],X.pred.orig = D.test$X, X01=T,y.var = T,verbose = F,NN = NN, parallel = parallel)
      rmse_flagp[i+1] = Rtools:::get_rmse(flagp.pred$y.mean,D.test$Y)
      mape_flagp[i+1] = Rtools:::get_mape(flagp.pred$y.mean,D.test$Y)
      force_refit = T
      cat('rmse new:',rmse_flagp[i+1],'\n')
    }
  }

  if(plot_during_seqd & i %% 5 == 0){
    par(mfrow=c(1,3))
    plot(m_init:(m_init+i),imse[1:(i+1)],type='l',col='cornflowerblue',ylim=c(0,max(c(imse,imse_flagp))),ylab='IMSE',xlab='M')
    lines(m_init:(m_init+i),imse_flagp[1:(i+1)],col='maroon')
    plot(m_init:(m_init+i),rmse[1:(i+1)],type='l',col='cornflowerblue',ylim=c(0,max(c(rmse,rmse_flagp))),ylab='RMSE',xlab='M')
    lines(m_init:(m_init+i),rmse_flagp[1:(i+1)],col='maroon')
    # plot(m_init:(m_init+i),mape[1:(i+1)],type='l',col='cornflowerblue',ylim=c(0,max(c(mape,mape_flagp))),ylab='MAPE',xlab='M')
    # lines(m_init:(m_init+i),mape_flagp[1:(i+1)],col='maroon')
    plot((m_init+1):(m_init+i),cumsum(lagp.time[1:i]),type='l',col='cornflowerblue',ylim=c(0,max(sum(flagp.time),sum(lagp.time))),ylab='Cumulative time (s)',xlab='M')
    lines((m_init+1):(m_init+i),cumsum(flagp.time[1:i]),col='maroon')
  }
  if(save_results)
    save(imse,imse_flagp,rmse,rmse_flagp,mape,mape_flagp,
         m_init,lagp.time,flagp.time,flagp.design.time,seed,
         file=save_file)
}

if(plot_during_seqd){
  par(mfrow=c(1,3),mar=c(5,5,4,4))
  plot((0:(length(imse)-1))+m_init,imse,type='l',col='cornflowerblue',ylim=c(0,max(c(imse,imse_flagp))),ylab='IMSE',xlab='M',lwd=2)
  lines((0:(length(imse)-1))+m_init,imse_flagp,col='maroon',lwd=2)
  legend('topright',inset=c(.05,.05),legend=c('laGP-alc','FlaGP-nn'),lwd=2,col=c('cornflowerblue','maroon'),lty=1,cex=1.5)
  # plot((0:(length(imse)-1))+m_init,rmse,type='l',col='cornflowerblue',ylim=c(0,.2),ylab='RMSE',xlab='M',lwd=2)
  # lines((0:(length(imse)-1))+m_init,rmse_flagp,col='maroon',lwd=2)
  plot((0:(length(imse)-1))+m_init,mape,type='l',col='cornflowerblue',ylim=c(0,2),ylab='MAPE',xlab='M',lwd=2)
  lines((0:(length(imse)-1))+m_init,mape_flagp,col='maroon',lwd=2)
  legend('topright',inset=c(.05,.05),legend=c('laGP-alc','FlaGP-nn'),lwd=2,col=c('cornflowerblue','maroon'),lty=1,cex=1.5)
  plot((0:(length(cumsum(lagp.time))-1))+m_init,log(cumsum(lagp.time),base = 10),type='l',lwd=2,col='cornflowerblue',
       ylim=c(log(flagp.time[1],base=10),max(log(sum(flagp.time),base=10),log(sum(lagp.time),base=10))),ylab='log Cumulative time (s)',xlab='M')
  lines((0:(length(cumsum(lagp.time))-1))+m_init,log(cumsum(flagp.time),base = 10),col='maroon',lwd=2)
  legend('topleft',inset=c(.05,.05),legend=c('laGP-alc','FlaGP-nn'),col=c('cornflowerblue','maroon'),lwd=2,lty=1,cex=1.5)
}


if(save_results)
  save(imse,imse_flagp,rmse,rmse_flagp,mape,mape_flagp,
       m_init,lagp.time,flagp.time,flagp.design.time,seed,
       file=save_file)

