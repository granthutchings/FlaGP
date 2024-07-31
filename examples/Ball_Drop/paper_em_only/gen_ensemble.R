t_at_d = function(C,d,R,g){
  (acosh(exp(C*d/R)))/(sqrt(C*g/R))
}

# number of x inputs
p.x = 2
# number of simulations
m = c(98,242,578,1058,5618,10082)

# physically relevant bounds on X
X.range = matrix(c(.1,2,.01,.5),nrow=p.x,ncol=2,byrow = T)
# distances where time is evaluated
y.ind.sim = as.matrix(seq(1,25,1)) # distance (meters)
n.y = length(y.ind.sim)

# sim data
set.seed(11)
reps = 100

X = list()
Y = list()


for(j in 1:length(m)){
  file = paste0('~/FlaGP/examples/Ball_Drop/paper_em_only/ensemble_m_',m[j],'.RData')
  M = m[j]
  for(k in 1:reps){
    # create design over (R,C) on [0,1]
    X[[k]] = lhs::create_oalhs(m[j],p.x,T,F)

    # put inputs on physically relavant scale
    X[[k]][,1] = X[[k]][,1] * (X.range[1,2]-X.range[1,1]) + X.range[1,1]
    X[[k]][,2] = X[[k]][,2] * (X.range[2,2]-X.range[2,1]) + X.range[2,1]
    # pairs(X[plot.sample,],labels = c('R','C'))

    # generate simulations
    Y[[k]] = matrix(nrow=n.y,ncol=m[j])
    for(i in 1:m[j]){
      Y[[k]][,i] = t_at_d(C=X[[k]][i,2],d=y.ind.sim,R=X[[k]][i,1],g=9.8)
    }
  }
  save(X,Y,y.ind.sim,n.y,M,file = file)
}

# generate prediction data on a lattice with boundaries at .05 and .95 quantiles
# of the training ensemble so we don't have to worry about extrapolation
R.pred.range = qunif(c(.05,.95),X.range[1,1],X.range[1,2])
C.pred.range = qunif(c(.05,.95),X.range[2,1],X.range[2,2])

X.pred = mined::Lattice(101,2)
n.pred = nrow(X.pred)
X.pred[,1] = X.pred[,1] * (R.pred.range[2]-R.pred.range[1]) + R.pred.range[1]
X.pred[,2] = X.pred[,2] * (C.pred.range[2]-C.pred.range[1]) + C.pred.range[1]
pairs(X.pred)
apply(X.pred,2,range)

Y.pred = matrix(nrow=n.y,ncol=n.pred)
for(i in 1:n.pred){
  Y.pred[,i] = t_at_d(C=X.pred[i,2],d=y.ind.sim,R=X.pred[i,1],g=9.8)
}
matplot(Y.pred,type='l')

save(X.pred,Y.pred,file='~/FlaGP/examples/Ball_Drop/paper_em_only/test_data.RData')

