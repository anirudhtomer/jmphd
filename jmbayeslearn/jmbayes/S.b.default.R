S.b.default <-
function (t, b, ii, Mats, log = FALSE) {
    if (t == 0) {
        if (log) return(0) else return(0)
    }

  id.i <- id %in% ii
  idT.i <- idT %in% ii
  ids.i <- ids %in% ii
  
  wk <- Mats$wk
  
  W_time_independent <- if(!is.null(W)){
    W[ii, , drop = FALSE]
  }else{
    rep(0, length(gammas.new))
  }
  
  ind <- Mats$st < min(Mats$kn)
  wk[ind] <- 0
  log.survival = -logSurvival_cpp(Mats$Xs, betas.new,  Mats$Zs, b, alphas.new, W_time_independent, 
                                  gammas.new, Mats$W2s, Bs.gammas.new, wk, Mats$P)
    

  if (log) log.survival else exp(log.survival)
}
