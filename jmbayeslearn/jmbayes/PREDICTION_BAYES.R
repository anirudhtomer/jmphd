survfitJM.JMbayes <- function (object, newdata, type = c("SurvProb", "Density"), 
                               idVar = "id", simulate = TRUE, survTimes = NULL, 
                               last.time = NULL, LeftTrunc_var = NULL, M = 200L, 
                               CI.levels = c(0.025, 0.975), 
                               log = FALSE, scale = 1.6, weight = rep(1, nrow(newdata)), 
                               init.b = NULL, seed = 1L, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0L)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    
    #Type tells you whether its survival or longitudinal prediction
    type <- match.arg(type)
    
    #TT tells the survival time for all the objects`
    TT <- object$y$Time

    # Check if the survival probabilities are to be calculated at times given by user
    if (is.null(survTimes) || !is.numeric(survTimes)) {
        survTimes <- seq(min(TT), quantile(TT, 0.90) + 0.01, length.out = 35L)
    }
    if (type != "SurvProb") simulate <- TRUE
    
    # which variable in the model corresponds to the time
    timeVar <- object$timeVar

    df.RE <- object$y$df.RE

    # td-both, td-value or td-extra
    param <- object$param

    #which density does the response variable follow
    densLong <- object$Funs$densLong

    # No clue what's this
    hasScale <- object$Funs$hasScale

    #Any left truncation in the data?
    anyLeftTrunc <- object$y$anyLeftTrunc

    # distribution for the random effects
    densRE <- object$Funs$densRE

    # the transformation functions in the cox model. x^2 and x*some other variable
    transFun.value <- object$Funs$transFun.value
    transFun.extra <- object$Funs$transFun.extra
    extraForm <- object$Forms$extraForm

    #don't care, some indexes 
    indFixed <- extraForm$indFixed
    indRandom <- extraForm$indRandom
    indBetas <- object$y$indBetas

    #longitudinal covariates X and Z
    TermsX <- object$Terms$termsYx
    TermsZ <- object$Terms$termsYz
    TermsX.extra <- object$Terms$termsYx.extra
    TermsZ.extra <- object$Terms$termsYz.extra

    #mfX has reponse and fixed effects, mfZ has random effects of the new data set of all the patients
    mfX <- model.frame.default(TermsX, data = newdata)
    mfZ <- model.frame.default(TermsZ, data = newdata)

    #formula Y to X, and Y to Z
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$Forms$formYz

    #the weight function if you use weighted integral for association function
    estimateWeightFun <- object$estimateWeightFun
    weightFun <- object$Funs$weightFun

    #max of survival time of all subjects
    max.time <- max(TT)


    # can't figure out what's happening here
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
        rep(TRUE, nrow(newdata))
    } else {
        !seq_len(nrow(newdata)) %in% na.ind
    }
    id <- as.numeric(unclass(newdata[[idVar]]))
    id <- id. <- match(id, unique(id))
    id <- id[na.ind]


    # response and covariates of the model we are interested in. Y is longtudinal
    y <- model.response(mfX)
    X <- model.matrix.default(formYx, mfX)
    Z <- model.matrix.default(formYz, mfZ)[na.ind, , drop = FALSE]


    #data.id is for time to event information of the patients for whom prediction is being done
    TermsT <- object$Terms$termsT
    data.id <- newdata[tapply(row.names(newdata), id, tail, n = 1L), ]
    
    #This one replicates each row of data.id as many times as gaussian quadrature requires
    data.s <- data.id[rep(1:nrow(data.id), each = object$control$GQsurv.k), ]
    
    # get the ID number of the subject into idT
    idT <- data.id[[idVar]]
    idT <- match(idT, unique(idT)) #here convert ID's into a a series of 1,2,3
    
    ids <- data.s[[idVar]]
    ids <- match(ids, unique(ids))

    if (type != "SurvProb") {
        SurvT <- model.response(model.frame(TermsT, data.id))
        Time <- SurvT[, 1]
        event <- SurvT[, 2]
    }

    #mfT means model and its covariates for survival part of the model
    mfT <- model.frame.default(delete.response(TermsT), data = data.id)
    formT <- if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else {
        tt <- attr(delete.response(TermsT), "term.labels")
        if (length(tt)) reformulate(tt) else reformulate("1")
    }

    #take the model mfT and formula formT and get the data matrix W
    W <- model.matrix.default(formT, mfT)[, -1L, drop = FALSE]
    
    #take the column of times where longitudinal measurements were taken and using split function 
    #divide them into separate vectors for every subject
    obs.times <- split(newdata[[timeVar]][na.ind], id)

    #the last longitudinal measurement time for every subject
    last.time <- if (is.null(last.time)) {
        tapply(newdata[[timeVar]], id., tail, n = 1L)
    } else if (is.character(last.time) && length(last.time) == 1L) {
        tapply(newdata[[last.time]], id., tail, n = 1L)
    } else if (is.numeric(last.time)) {
        rep_len(last.time, length.out = nrow(data.id))
    } else {
        stop("\nnot appropriate value for 'last.time' argument.")
    }

    #all times that are greater than the last longitudinal measurement time 
    #because if he lived till last measurement then for sure he is alive
    times.to.pred <- if (type == "SurvProb") {
        lapply(last.time, function (t) survTimes[survTimes > t])
    } else {
        as.list(Time)
    }


    TimeL <- if (!is.null(anyLeftTrunc) && anyLeftTrunc) {
        if (is.null(LeftTrunc_var) || is.null(newdata[[LeftTrunc_var]])) {
            warning("The original joint model was fitted in a data set with left-",
                    "truncation and\nargument 'LeftTrunc_var' of survfitJM() has not ", 
                    "been specified.\n")
        }
        TimeL <- newdata[[LeftTrunc_var]]
        tapply(TimeL, id, head, n = 1)
    }

    #Get the means of all of the parameters
    n <- length(TT)
    #n.tp corresponds to the number of subjects
    n.tp <- length(last.time)

    #number of columns in X, Z and W matrices of the joint model
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)


    if (ncww == 0L)
        W <- NULL
    lag <- object$y$lag
    betas <- object$postMeans$betas
    sigma <- object$postMeans$sigma
    D <- object$postMeans$D
    gammas <- object$postMeans$gammas
    alphas <- object$postMeans$alphas
    Dalphas <- object$postMeans$Dalphas
    shapes <- object$postMeans$shapes
    Bs.gammas <- object$postMeans$Bs.gammas
    list.thetas <- list(betas = betas, sigma = sigma, gammas = gammas, alphas = alphas, 
                        Dalphas = Dalphas, shapes = shapes, Bs.gammas = Bs.gammas, D = D)
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    environment(log.posterior.b) <- environment(S.b) <- environment(logh.b) <- environment()
    environment(hMats) <- environment(ModelMats) <- environment()


    # The last longitudinal measurement time for each subject, in a list
    obs.times.surv <- split(data.id[[timeVar]], idT)
    
    #just created empty list...to be filled in the next for loop
    survMats <- survMats.last <- vector("list", n.tp) 

    #Run as many times as the number of subjects
    for (i in seq_len(n.tp)) {
        #lapply on every single time at which we want to make a prediction
        survMats[[i]] <- lapply(times.to.pred[[i]], ModelMats, ii = i, timeL = TimeL[i])
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i, timeL = TimeL[i])
    }
    if (type != "SurvProb")
        hazMats <- hMats(Time)

    # calculate the Empirical Bayes estimates and their (scaled) variance
    modes.b <- matrix(0, n.tp, ncz)
    invVars.b <- Vars.b <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        betas.new <- betas
        sigma.new <- sigma
        D.new <- D
        gammas.new <- gammas
        alphas.new <- alphas
        Dalphas.new <- Dalphas
        shapes.new <- shapes
        Bs.gammas.new <- Bs.gammas
        ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y, Mats = tt, ii = i)
        start <- if (is.null(init.b)) rep(0, ncz) else init.b[i, ]
        opt <- try(optim(start, ff, y = y, tt = survMats.last, i = i, 
            method = "BFGS", hessian = TRUE), silent = TRUE)
        
        #if error happened during optimization routine
        if (inherits(opt, "try-error")) {
            gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, i = i)
            opt <- optim(start, ff, gg, y = y, tt = survMats.last, 
                i = i, method = "BFGS", hessian = TRUE, 
                control = list(parscale = rep(0.1, ncz)))
        } 
        modes.b[i, ] <- opt$par
        invVars.b[[i]] <- opt$hessian/scale
        Vars.b[[i]] <- scale * solve(opt$hessian)
    }


    if (!simulate) {
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            S.last <- S.b(last.time[i], modes.b[i, ], i, survMats.last[[i]])
            S.pred <- numeric(length(times.to.pred[[i]]))
            for (l in seq_along(S.pred))
                S.pred[l] <- S.b(times.to.pred[[i]][l], modes.b[i, ], i, survMats[[i]][[l]])
            res[[i]] <- cbind(times = times.to.pred[[i]], predSurv = weight[i] * S.pred / S.last)
            rownames(res[[i]]) <- seq_along(S.pred) 
        }
    } else {
        set.seed(seed)
        out <- vector("list", M)
        success.rate <- matrix(FALSE, M, n.tp)
        b.old <- b.new <- modes.b
        if (n.tp == 1)
            dim(b.old) <- dim(b.new) <- c(1L, ncz)
        mcmc <- object$mcmc
        mcmc <- mcmc[names(mcmc) != "b"]
        if (M > nrow(mcmc$betas)) {
            warning("'M' cannot be set greater than ", nrow(mcmc$betas))
            M <- nrow(mcmc$betas)
            out <- vector("list", M)
            success.rate <- matrix(FALSE, M, n.tp)
        }
        samples <- sample(nrow(mcmc$betas), M)
        mcmc[] <- lapply(mcmc, function (x) x[samples, , drop = FALSE])
        proposed.b <- mapply(rmvt, mu = split(modes.b, row(modes.b)), Sigma = Vars.b, 
                             MoreArgs = list(n = M, df = 4), SIMPLIFY = FALSE)
        proposed.b[] <- lapply(proposed.b, function (x) if (is.matrix(x)) x else rbind(x))
        dmvt.proposed <- mapply(dmvt, x = proposed.b, mu = split(modes.b, row(modes.b)),
                                Sigma = Vars.b, MoreArgs = list(df = 4, log = TRUE), 
                                SIMPLIFY = FALSE)
        for (m in 1:M) {
            # Step 1: extract parameter values
            betas.new <- mcmc$betas[m, ]
            if (hasScale)
                sigma.new <- mcmc$sigma[m, ]
            if (!is.null(W))
                gammas.new <- mcmc$gammas[m, ]
            if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
                alphas.new <- mcmc$alpha[m, ]
            if (param %in% c("td-extra", "td-both"))
                Dalphas.new <- mcmc$Dalphas[m, ]
            if (estimateWeightFun)
                shapes.new <- mcmc$shapes[m, ]
            D.new <- mcmc$D[m, ]; dim(D.new) <- dim(D)
            Bs.gammas.new <- mcmc$Bs.gammas[m, ]
            if (type != "SurvProb") {
                logHaz <- logh.b(modes.b, hazMats)
            }
            SS <- vector("list", n.tp)
            for (i in seq_len(n.tp)) {
                # Step 2: simulate new random effects values
                p.b <- proposed.b[[i]][m, ]
                dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], invSigma = invVars.b[[i]], 
                                 df = 4, log = TRUE)
                dmvt.prop <- dmvt.proposed[[i]][m]
                a <- min(exp(log.posterior.b(p.b, y, survMats.last, ii = i) + dmvt.old - 
                        log.posterior.b(b.old[i, ], y, survMats.last, ii = i) - dmvt.prop), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (!is.na(ind) && ind)
                    b.new[i, ] <- p.b
                # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
                logS.last <- S.b(last.time[i], b.new[i, ], i, survMats.last[[i]], 
                                 log = TRUE)
                logS.pred <- numeric(length(times.to.pred[[i]]))
                for (l in seq_along(logS.pred))
                    logS.pred[l] <- S.b(times.to.pred[[i]][l], b.new[i, ], i, 
                                        survMats[[i]][[l]], log = TRUE)
                if (type != "SurvProb") {
                    logS.pred <- event[i] * logHaz[i] + logS.pred
                }
                SS[[i]] <- weight[i] * if (log) logS.pred - logS.last else exp(logS.pred - logS.last)
            }
            b.old <- b.new
            out[[m]] <- SS
        }
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            rr <- sapply(out, "[[", i)
            if (!is.matrix(rr))
                rr <- rbind(rr)
            res[[i]] <- cbind(
                times = times.to.pred[[i]],
                "Mean" = rowMeans(rr, na.rm = TRUE),
                "Median" = apply(rr, 1L, median, na.rm = TRUE),
                "Lower" = apply(rr, 1L, quantile, probs = CI.levels[1], na.rm = TRUE),
                "Upper" = apply(rr, 1L, quantile, probs = CI.levels[2], na.rm = TRUE)
            )
            rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
        }
    }
    y <- split(y, id)
    newdata. <- do.call(rbind, mapply(function (d, t) {
        d. <- rbind(d, d[nrow(d), ])
        d.[[timeVar]][nrow(d.)] <- t
        d.
    }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
    id. <- as.numeric(unclass(newdata.[[idVar]]))
    id. <- match(id., unique(id.))
    mfX. <- model.frame(delete.response(TermsX), data = newdata.)
    mfZ. <- model.frame(TermsZ, data = newdata.)
    X. <- model.matrix(formYx, mfX.)
    Z. <- model.matrix(formYz, mfZ.)
    fitted.y <- split(c(X. %*% betas) + rowSums(Z. * modes.b[id., , drop = FALSE]), id.)
    names(res) <- names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
    res <- list(summaries = res, survTimes = survTimes, last.time = last.time, 
        obs.times = obs.times, y = y, 
        fitted.times = split(newdata.[[timeVar]], factor(newdata.[[idVar]])), 
        fitted.y = fitted.y, ry = range(object$y$y, na.rm = TRUE),
        nameY = paste(object$Forms$formYx)[2L], modes.b = modes.b)
    if (simulate) {
        res$full.results <- out
        res$success.rate <- success.rate
    }
    if (simulate) rm(list = ".Random.seed", envir = globalenv())
    class(res) <- "survfit.JMbayes"
    res
}
