# Using Bayes rule
# f(b|T>t, Y(t), theta) = f(T>t, Y(t) | b, theta) f(b|theta) f(theta) / f(T>t, Y(y))
# However in empirical bayes you won't consider f(theta) but will rather use MLE of theta
# Thus f(b|T>t, Y(t), theta) = f(T>t| b, theta) f(Y(t)| b, theta) f(b|MLE(theta)) / f(T>t, Y(y))

log.posterior.b.normal <- function (b, y, Mats, ii) {
  id.i <- id %in% ii
  idT.i <- idT %in% ii
  ids.i <- ids %in% ii
  X.i <- X[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  
  log.p.yb <- dnorm_cpp(y[id.i], X.i, betas.new, Z.i, b, sigma.new)
  log.p.b <- dmvnorm_cpp(t(b), t(rep(0, ncol(Z.i))), D.new, logd = T)
  
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
  
  tt <- if(estimateWeightFun){
    wFun <- wk2 * weightFun(st2, shapes.new, max.time)
    Yu <- transFun.value(P2 * fastSumID(wFun * c(Xu %*% betas.new + Zu %*% b), id.GK2),
                         data.s[ids.i, ])
    c(as.matrix(Yu) %*% alphas.new)
  }else{
    if(param=="td-both"){
      Ys <- transFun.value(XbetaZb_cpp(Xs, betas.new,Zs, b), data.s[ids.i, ])
      Ys.extra <- transFun.extra(XbetaZb_cpp(Xs.extra, betas.new[indFixed],
                                         Zs.extra, b[indRandom]), data.s[ids.i, ])
      XbetaZb_cpp(as.matrix(Ys), alphas.new, 
              as.matrix(Ys.extra), Dalphas.new)
    }else if(param=="td-value"){
      Ys <- transFun.value(XbetaZb_cpp(Xs, betas.new,Zs, b), data.s[ids.i, ])
      Xbeta_cpp(as.matrix(Ys), alphas.new)
    }else if(param=="td-extra"){
      Ys.extra <- transFun.extra(XbetaZb_cpp(Xs.extra, betas.new[indFixed],
                                         Zs.extra, b[indRandom]), data.s[ids.i, ])  
      Xbeta_cpp(as.matrix(Ys.extra), Dalphas.new)
    }else if(param=="shared-betasRE"){
      rep(sum((betas[indBetas] + b) * alphas.new), length(st))
    }else if(param=="shared-RE"){
      rep(sum(b * alphas.new), length(st))
    }
  }
  
  #the part of W matrix which doesn't depend on time
  eta.tw = if (!is.null(W)) {
    as.vector(Xbeta_cpp(W[ii, , drop = FALSE], gammas.new))
  } else {
    0
  }
  
  #log of baseline hazard or log(h0(t)) = W2s %*% Bs.gammas.new. The "Bs.gammas.new" are the coefficients in spline approx
  Vi <- exp(c(Xbeta_cpp(W2s, Bs.gammas.new)) + tt)
  log.survival <- - sum(exp(eta.tw) * P * fastSumID(wk * Vi, idT))
  if (all(st == 0))
    log.survival <- 1
  
  log.p.yb + log.survival + log.p.b
}
