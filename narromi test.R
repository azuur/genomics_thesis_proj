library(rbenchmark)

n=200; d=50

res <- benchmark("LADLasso" = {
  X=matrix(rnorm(n*d), nrow=n, ncol=d)
  beta=c(rep(2,6), rep(0, 44))
  y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
  output.LADLasso=LADlasso(y,X, beta.ini=LAD(y, X), lambda=0.2, adaptive=F)
  beta.est=output.LADLasso$beta
},
"LP_TGN" = {
  X=matrix(rnorm(n*d), nrow=n, ncol=d)
  beta=c(rep(2,6), rep(0, 44))
  y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
  b <- LP_TGN(y,X,0.2*nrow(X))
},
replications = 1000,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))

res