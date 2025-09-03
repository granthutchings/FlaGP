library(FlaGP)
library(Rtools)
library(doParallel)
load = FALSE
load_data = TRUE
loadfile = 'examples/Aluminum/seq_design/imse_8_PC_refit_50_test_set_no_boot_results.RData'
if(load){
  load(loadfile)
  start = i + 1
} else{
  if(load_data){
    load('examples/Aluminum/seq_design/flyer_data_seq_design.RData')
  } else{
    inputs = read.csv('examples/Aluminum/data/flyer_sim_inputs.csv')
    input_names = c('A', 'B', 'C', 'n', 'm', 'v1', 'v2', 'v3', 'G1', 'del2', 'del3')

    sims104 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/104_sim_aligned.csv',header = F))
    sims105 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/105_sim_aligned.csv',header = F))
    sims106 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/106_sim_aligned.csv',header = F))
    outputs = cbind(sims104,sims105,sims106)
    rm(sims104,sims105,sims106)
    obs104 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/104_obs_aligned.csv',header = F))
    obs105 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/105_obs_aligned.csv',header = F))
    obs106 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/106_obs_aligned.csv',header = F))
    obs_all = matrix(c(obs104,obs105,obs106))
    rm(obs104,obs105,obs106)
    time104 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/104_time.csv',header = F))
    time105 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/105_time.csv',header = F))
    time106 = as.matrix(read.csv('examples/Aluminum/data/python_clean_data/106_time.csv',header = F))
    time_all = time104
    time_all = c(time_all,time105-min(time105)+max(time_all))
    time_all = c(time_all,time106-min(time106)+max(time_all))
    rm(time104,time105,time106)

    M = nrow(inputs)
    plot.id = as.integer(seq(1,M,length.out=100))
    matplot(time_all,t(outputs[plot.id,]),type='l',col='darkorange')
    lines(time_all,obs_all,lwd=2,col='black')

    # start with random subset
    set.seed(4321)
    sub = sample(1:nrow(inputs),2000+110+2000+2000,replace = F)
    inputs = inputs[sub,]
    inputs = FlaGP:::transform_xt(inputs)$sim$X$trans
    outputs = outputs[sub,]

    # create an initial set of points using maximin space filling
    dt = ncol(inputs)
    m_init = 10*dt
    X_init = Rtools:::maximin_from_X_2(m_init,inputs,scale = F)
    # pairs(X_init$lhs)
    Y_init = outputs[X_init$orig_index,]
    outputs = outputs[-X_init$orig_index,]

    # now a set of integration points
    m_int = 2000
    X_int = Rtools:::maximin_from_X_2(m_int,X_init$newX,scale = F)
    pairs(X_int$lhs)
    Y_int = outputs[X_int$orig_index,]
    outputs = outputs[-X_int$orig_index,]

    # a candidate set
    m_cand = 2000
    X_cand = Rtools:::maximin_from_X_2(m_cand,X_int$newX,scale = F)
    pairs(X_cand$lhs)
    Y_cand = outputs[X_cand$orig_index,]
    outputs = outputs[-X_cand$orig_index,]

    # a test set
    m_test = 2000
    X_test = X_cand$newX
    pairs(X_test)
    Y_test = outputs
    # outputs = outputs[-X_test$orig_index,]
    rm(outputs)
  }
  ### emulation ###

  # define model and fit emulator
  model.em = flagp(X.sim = X_init$lhs, Y.sim = t(Y_init), y.ind.sim = matrix(time_all),
                  transform_x = F,
                  pct.var = .99, seed = 123,
                  ls.subsample = 'strat', ls.subsample.size = 1024, ls.K = 1,
                  nug.est = T, ls.parallel = T, rsvd = F, ls.nugget = 1e-8) # not estimating a nugget can cause problems for the fitting procedure
  model.em$lengthscales$XT[[1]]
  cumsum(model.em$basis$sim$pct.var)

  # this makes me wonder if we should take the active learning approach to basis vectors. SVD is not great here,
  # the responses clearly violate the normal assumption
  # plot(model.calib,sim.subsample = sample(1:10000,1000))
  plot(model.em)

  ### Test set results - emulator only prediction
  NN = 50
  int.pred = predict(model.em,X.pred.orig = X_test,end.eta = NN,y.conf.int = T)
  int.rmse = Rtools:::get_rmse(int.pred$y.mean,t(Y_test))

  # sequential design
  n_points_to_add = 500
  refitevery = 50
  plotevery = 5
  n.pc.sd = 1
  method = 'maxvar'
  flagp.model = list()
  flagp.model[[1]] = model.em
  xcandflagp = X_cand$lhs
  xtint = X_int$lhs
  time = imse = rmse = numeric(n_points_to_add+1)
  rmse[1] = int.rmse
  imse[1] = mean(int.pred$y.var)
  start = 1
}

parallel = T
if(parallel){
  cl = parallel::makeCluster(model.em$basis$sim$n.pc)
  doParallel::registerDoParallel(cl)
}

last_save_time = Sys.time()
cat('\nsequential design using method:',method,'\n')
for(i in start:n_points_to_add){
  cat(i,' ')
  m = m_init + i
  refit = ifelse(m%%refitevery==0,T,F)
  ####################################################################################################
  # FlaGP
  ptm = proc.time()
  if(i==1)
    XintVars = NULL

  # find new point
  xnewF = FlaGP::seq_design(model = flagp.model[[max(1,i-1)]],
                            Xcand = xcandflagp,
                            Xint = xtint,
                            XintPredVars = XintVars,
                            end = NN,method=method,n.pc = n.pc.sd,
                            verbose=T, parallel = parallel)


  # remove the row from xcand and XcandVars
  XintVars = xnewF$XintPredVars
  xcandflagp = xcandflagp[-xnewF$cand_id,]

  # find y associated with new point
  idnew = Rtools:::row_match(matrix(xnewF$X_new,nrow=1),X_cand$lhs)
  yneworig = Y_cand[idnew$rowid,]

  # update model with new point
  k = 1
  flagp.model[[i]] = FlaGP::flagp_update(flagp.model[[max(1,i-1)]],matrix(xnewF$X_new,nrow=1),Ynew_orig = yneworig,refit = refit,
                                         seed=k*i,n.pc=flagp.model[[max(1,i-1)]]$basis$sim$n.pc,parallel=parallel,make.cluster=F)

  # flag candidate points that need variance recalculated
  if(refit | method == 'maxvar'){
    # the model was refit so we need to recalculate everything
    XintVars = NULL
  } else if(method == 'imse'){
    XintVars[xnewF$effected] = NA
  }
  time[i+1] = -(ptm - proc.time())[3]

  flagp.pred = predict(flagp.model[[i]],X.pred.orig = X_test, y.var = T,verbose = F,end.eta = NN, parallel = T)
  imse[i+1] = mean(flagp.pred$y.var)
  rmse[i+1] = Rtools:::get_rmse(flagp.pred$y.mean,t(Y_test))
  cat('rmse prev:',rmse[i],' rmse new:',rmse[i+1],'\n')
  
  # if(refit & rmse[i+1]>rmse[i]){
  #   while(rmse[i+1]>rmse[i]){
  #     k = k + 1
  #     cat('refitting with seed',k*i,'\n')
  #     flagp.model[[i]] = FlaGP::flagp_update(flagp.model[[max(1,i-1)]],matrix(xnewF$X_new,nrow=1),Ynew_orig = yneworig,refit = refit,
  #                                            seed=as.integer(k*i),n.pc=flagp.model[[max(1,i-1)]]$basis$sim$n.pc,parallel=parallel,make.cluster=F)
  #     flagp.pred = predict(flagp.model[[i]],X.pred.orig = X_test, y.var = T,verbose = F,end.eta = NN, parallel = T)
  #     imse[i+1] = mean(flagp.pred$y.var)
  #     rmse[i+1] = Rtools:::get_rmse(flagp.pred$y.mean,t(Y_test))
  #     cat('rmse prev:',rmse[i],' rmse new:',rmse[i+1],'\n')
  #   }
  # }
  if(i%%plotevery==0){
    par(mfrow=c(1,2),mar=c(4,4,4,4))
    plot(m_init:(m_init+i),rmse[1:(i+1)],type='l',col='dodgerblue',lwd=2,xlab='M',ylab='RMSE')
    plot(m_init:(m_init+i),imse[1:(i+1)],type='l',col='darkorange',lwd=2,xlab='M',ylab='IMSPE')
  }
  current_time <- Sys.time()
  if(m%%50==0 | i==n_points_to_add | difftime(current_time, last_save_time, units = "secs") >= 3600){
    last_save_time = current_time
    cat('saving current status \n')
    if(method %in% c('maximin','maximin-scaled')){
      save.image(file = paste0('examples/Aluminum/seq_design/',method,'_refit_',refitevery,'_test_set_no_boot_results_nugget.RData'))
    } else{
      save.image(file = paste0('examples/Aluminum/seq_design/',method,'_',n.pc.sd,'_PC_refit_',refitevery,'_test_set_no_boot_results_nugget.RData'))
    }
  }
}
