dat=data.frame(x=c(1,2,3,4,5,6), 
               y=c(1,3,5,6,8,12))

min.RSS <- function(par,data, Mats, ii, b, y) {
  
  id.i <- id %in% ii
  idT.i <- idT %in% ii
  ids.i <- ids %in% ii
  X.i <- X[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  mu.y <- as.vector(X.i %*% betas.new + Z.i %*% b)
  logY <- densLong(y[id.i], mu.y, sigma.new, log = TRUE, data = newdata[id.i, ])
  log.p.yb <- sum(logY)
  log.p.b <- densRE(b, mu = rep(0, ncol(Z.i)), D = D.new, log = TRUE, prop = FALSE)
  MM <- Mats[[ii]]
  st <- MM$st
  st2 <- MM$st2
  wk <- MM$wk
  wk2 <- MM$wk2
  P <- MM$P
  P2 <- MM$P2
  W2s <- MM$W2s
  Xs <- MM$Xs
  Zs <- MM$Zs
  Xs.extra <- MM$Xs.extra
  Zs.extra <- MM$Zs.extra
  Xu <- MM$Xu
  Zu <- MM$Zu
  ind <- MM$ind
  idT <- MM$idT
  id.GK2 <- MM$id.GK2
  
  
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
  matrix(rnorm(20*6, sd=10), nrow=20, ncol=6) %*% rnorm(6, sd = 100) 
    
    with(data, sum((par[1] + par[2] * x - y)^2))
}

C = microbenchmark({
ret = foreach(i=1:5, .export = c("dmvnorm")) %dopar%{
  opt = try(optim(par = c(0, 1), min.RSS, data = dat, hessian = T))
  
  x = 5
  
  if(inherits(opt, "try-error")){
    x = 10
  }
  
  x
}})