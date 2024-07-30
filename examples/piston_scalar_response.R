# This example illustrates how to define a scalar response model for FlaGP

#### We start with a simple emulation problem on the piston function
library(duqling)
set.seed(1)
X = tgp::lhs(1000,matrix(rep(c(0,1),7),nrow=7,byrow=T))
Y = apply(X,1,piston,scale01=T)
set.seed(11)
X.pred = tgp::lhs(1000,matrix(rep(c(0,1),7),nrow=7,byrow=T))

# rmse can be improved significantly at the cost of compute time
# by increasing ls.subsample.size
data = flagp(X.sim = X,Y.sim = Y,
             ls.subsample = 'strat', ls.subsample.size = 100,
             ls.m = 4, ls.K = 10, ls.nugget = 1e-8)

flagp.pred = predict(data,end.eta = 100,
                     X.pred.orig=X.pred,
                     n.samples = 1000,
                     conf.int = T)
Y.pred = apply(X.pred,1,piston,scale01=T)
rmse = sqrt(mean((Y.pred - flagp.pred$y.mean)^2))

### Calibration
data$lengthscales$XT # 2nd input is important, lets calibrate it
XT.obs = tgp::lhs(25,matrix(rep(c(0,1),7),nrow=7,byrow=T))
XT.obs[,2] = .25 # all obs data generated with x2=.5
Y.obs = apply(XT.obs,1,piston,scale01=T)
Y.obs = Y.obs + rnorm(length(Y.obs),0,.01)

data = flagp(X.sim = X[,-2],T.sim = X[,2,drop=F],Y.sim=Y,
             X.obs = XT.obs[,-2],Y.obs = Y.obs,
             ls.subsample = 'strat', ls.subsample.size = 100,
             ls.m = 4, ls.K = 10, ls.nugget = 1e-8)
calib = FlaGP::mcmc(data,n.samples = 10000,n.burn = 0,t.init = .25)
plot(calib)
