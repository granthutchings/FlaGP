load('~/FlaGP/examples/Ball_Drop/paper_em_only/test_data.RData')

# constants
m =     c(98,242,578,1058,5618,10082)
d_est = c(2,   4,  8,  16,  50,  100)
r_est = c(25, 25, 25,  25,   1,    1)
nugget = 1e-8
n.pred = nrow(X.pred)
reps = 100
n.y = 25

# results
interval.score = array(NA,dim=c(length(m),reps,n.y,n.pred,3))
MAPE = array(NA,dim=c(length(m),reps,n.y,n.pred))

for(j in 1:1){
  cat(m[j],'\n')
  file = paste0('~/FlaGP/examples/Ball_Drop/paper_em_only/ensemble_m_',m[j],'.RData')
  load(file)
  flagp.model = flagp.pred = list()
  for(i in 1:reps){
    if(i%%10==0){cat(i,' ')}
    flagp.model[[i]] = flagp(X.sim = X[[i]], Y.sim = Y[[i]], n.pc = 2,
                             ls.subsample = 'blhs', ls.m = d_est[j],
                             verbose=F, ls.nugget = 1e-8, ls.K = r_est[j])
    flagp.pred[[i]] = predict(flagp.model[[i]],X.pred.orig = X.pred, verbose = F,
                              n.samples = 100, conf.int = T)
    IS = FlaGP:::interval_score(Y.pred,y.conf.int = flagp.pred[[i]]$y.conf.int, terms = T)
    interval.score[j,i,,,1] = IS$i_score
    interval.score[j,i,,,2] = IS$width
    interval.score[j,i,,,3] = IS$penalty
    MAPE[j,i,,] = 100 * abs(flagp.pred[[i]]$y.mean - Y.pred) / Y.pred
  }
  file = paste0('~/FlaGP/examples/Ball_Drop/paper_em_only/flagp_results_m_',m[j],'.RData')
  IS = interval.score[j,,,,]
  mape = MAPE[j,,,]
  save(flagp.model,flagp.pred,IS,mape,file=file)
}

mean(MAPE,na.rm = T)
median(MAPE,na.rm = T)
mean(interval.score[,,,,1],na.rm = T)
median(interval.score[,,,,1],na.rm = T)
