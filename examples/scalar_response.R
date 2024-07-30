library(duqling)
set.seed(1)
X = tgp::lhs(10000,matrix(rep(c(0,1),7),nrow=7,byrow=T))
Y = apply(X,1,piston,scale01=T)
set.seed(11)
X.pred = tgp::lhs(1000,matrix(rep(c(0,1),7),nrow=7,byrow=T))

library(FlaGP)
fit = GpGp::fit_model(Y,X,covfun_name = 'matern35_scaledim')
COV = GpGp::matern35_isotropic(fit$covparms,X)
data = flagp(X,NULL,NULL,NULL,Y,
             ls.subsample.size = 500, ls.m = 4, ls.K = 10, ls.nugget = 1e-8)

flagp.pred = predict(data,end.eta = 100,
                     X.pred.orig=X.pred,
                     n.samples = 1000,
                     conf.int = T)
Y.pred =apply(X.pred,1,piston,scale01=T)
rmse = sqrt(mean((Y.pred - flagp.pred$y.mean)^2))

### Scaled Vecchia
source('~/Desktop/Computer Experiments papers/Reading Group/scaledVecchia/vecchia_scaled.R')
fit.sv=fit_scaled(Y,X,nug = 1e-8,nu=3.5,print.level = 1)
preds=predictions_scaled(fit=fit.sv,locs_pred=X.pred,m=100)
rmse.sc = sqrt(mean((Y.pred - preds)^2))
