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
  
  Ys <- transFun.value(MM$Xs %*% betas.new + MM$Zs %*% b)
  tt <- as.matrix(Ys) %*% alphas.new
  
  eta.tw <- if (!is.null(W)) {
    as.vector(W[ii, , drop = FALSE] %*% gammas.new)
  } else 0
  
  Vi <- exp(c(MM$W2s %*% Bs.gammas.new) + tt)
  log.survival <- - sum(exp(eta.tw) * MM$P * fastSumID(MM$wk * Vi, MM$idT))
  
  if (all(MM$st == 0))
    log.survival <- 1
  
  glob$logResList$y = c(glob$logResList$y, log.p.yb)
  glob$logResList$b = c(glob$logResList$b, log.p.b)
  glob$logResList$surv = c(glob$logResList$surv, log.survival)
  
  log.p.yb + log.survival + log.p.b
}
