S.b.allTimePoints <-
  function (t, b, ii, Mats, log = FALSE) {
    S.b.result =  vector("numeric", numberofPredictions)
    
    # if (t == 0) {
    #   if (log) return(0) else return(0)
    # }
    
    indTGreaterThan0 = c(1:numberofPredictions)[t!=0]
    
    ids.i <- ids %in% ii
    
    st <- Mats[[1]]$st
    
    #P <- Mats$P
    P.all <- sapply(Mats, FUN=function(x){x$P})
    
    P2 <- Mats$P2
    
    #Xs <- Mats$Xs
    Xs.all <- do.call(rbind, lapply(Mats, FUN=function(x){x$Xs}))
    #Zs <- Mats$Zs
    Zs.all <- do.call(rbind, lapply(Mats, FUN=function(x){x$Zs}))
    
    #Xs.extra <- Mats$Xs.extra
    Xs.extra.all <- do.call(rbind, lapply(Mats, FUN=function(x){x$Xs.extra}))
    
    #Zs.extra <- Mats$Zs.extra
    Zs.extra.all <- do.call(rbind, lapply(Mats, FUN=function(x){x$Zs.extra}))
    
    
    #W2s <- Mats$W2s
    W2s.all <- do.call(rbind, lapply(Mats, FUN=function(x){x$W2s}))
    
    tt.all <- if(estimateWeightFun){
      
      ######### TO BE COMPLETED ###########
      st2 <- Mats$st2
      wk2 <- Mats$wk2
      id.GK2 <- Mats$id.GK2
      Xu <- Mats$Xu
      Zu <- Mats$Zu
      
      wFun <- wk2 * weightFun(st2, shapes.new, max.time)
      Yu.all <- transFun.value(P2 * fastSumID(wFun * c(Xu %*% betas.new + Zu %*% b), id.GK2), 
                           data.s[ids.i, ])
      Xbeta_cpp(as.matrix(Yu.all), alphas.new)
      ######### TO BE COMPLETED ###########
    }else{
      data.s.all <- do.call(rbind, lapply(1:length(indTGreaterThan0), FUN=function(x){data.s[ids.i,]}))
      
      if(param=="td-both"){
        Ys.all <- transFun.value(XbetaZb_cpp(Xs.all, betas.new, Zs.all, b), data.s.all)
        Ys.extra.all <- transFun.extra(XbetaZb_cpp(Xs.extra.all, betas.new[indFixed], 
                                               Zs.extra.all, b[indRandom]), data.s.all)
        
        XbetaZb_cpp(as.matrix(Ys.all), alphas.new, 
                    as.matrix(Ys.extra.all), Dalphas.new)
      }else if(param=="td-value"){
        Ys.all <- transFun.value(XbetaZb_cpp(Xs.all, betas.new, Zs.all, b), data.s.all)
        Xbeta_cpp(as.matrix(Ys.all), alphas.new)
      }else if(param=="td-extra"){
        Ys.extra.all <- transFun.extra(XbetaZb_cpp(Xs.extra.all, betas.new[indFixed], 
                                                   Zs.extra.all, b[indRandom]), data.s.all)
        
        Xbeta_cpp(as.matrix(Ys.extra.all), Dalphas.new)
      }else if(param=="shared-betasRE"){
        rep(sum((betas[indBetas] + b) * alphas.new), length(st) * length(indTGreaterThan0))
      }else if(param=="shared-RE"){
        rep(sum(b * alphas.new), length(st) * length(indTGreaterThan0))
      }
    }
    
    eta.tw <- if (!is.null(W)) {
      as.vector(Xbeta_cpp(W[ii, , drop = FALSE], gammas.new))
    } else {
      0
    }
    eta.tw.all<-rep(eta.tw, length(indTGreaterThan0))
    
    Vi.all <- exp(Xbeta_cpp(W2s.all, Bs.gammas.new) + tt.all)
    
    fastSumID.res.all = vector("numeric", numberofPredictions)
    
    w2sLength = nrow(W2s.all)/numberofPredictions
    for(predictionNum in 1:numberofPredictions){
      ind <- Mats[[predictionNum]]$st < min(Mats[[predictionNum]]$kn)
      wk <- Mats[[predictionNum]]$wk
      wk[ind] <- 0
      Vi <- Vi.all[((predictionNum-1)*w2sLength + 1):(predictionNum*w2sLength)]
      
      fastSumID.res.all[predictionNum] <- fastSumID(wk * Vi, Mats[[predictionNum]]$idT)
    }
    
    S.b.result[indTGreaterThan0] <- -(exp(eta.tw.all) * P.all * fastSumID.res.all)
    
    # if (log){
    #   log.survival
    # }  else {
    #   exp(log.survival)
    # }
    
    if(log==F){
      S.b.result[indTGreaterThan0] = exp(S.b.result[indTGreaterThan0])
    }
    
    return(S.b.result)
  }
