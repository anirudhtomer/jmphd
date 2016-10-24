#no trans function, normal dist response, mv normal RE

log.posterior.b.default <- function (b, y, Mats, ii) {
  id.i <- id %in% ii
  idT.i <- idT %in% ii
  ids.i <- ids %in% ii
  X.i <- X[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  
  log.p.yb <- dnorm_cpp(y[id.i], X.i, betas.new, Z.i, b, sigma.new)
  log.p.b <- dmvnorm_cpp(t(b), t(rep(0, ncol(Z.i))), D.new, logd = T)
  
  
  MM <- Mats[[ii]]
  
  W_time_independent <- if(!is.null(W)){
    W[ii, , drop = FALSE]
  }else{
    rep(0, length(gammas.new))
  }
  
  log.survival = -logSurvival_cpp(MM$Xs, betas.new,  MM$Zs, b, alphas.new, W_time_independent, 
                              gammas.new, MM$W2s, Bs.gammas.new, MM$wk, MM$P)
                              
  if (all(MM$st == 0))
    log.survival <- 1
  
  log.p.yb + log.survival + log.p.b
}
