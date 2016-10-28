library("JMbayes")
library("lattice")
library("ggmcmc")
library(RcppArmadillo)
library(Rcpp)
library(inline)
library(profvis)
library(splines)
library(microbenchmark)
library(doParallel)

# cl<-makeCluster(detectCores())
# registerDoParallel(cl)
registerDoParallel(cores=detectCores())

prothro$t0 <- as.numeric(prothro$time == 0)
lmeFit.pro <- lme(pro ~ treat * (ns(time, 3) + t0),
                  random = list(id = pdDiag(form = ~ ns(time, 3))),
                  data = prothro)

coxFit.pro <- coxph(Surv(Time, death) ~ treat, data = prothros,
                    x = TRUE)

jointFit.pro <- jointModelBayes(lmeFit.pro, coxFit.pro, timeVar = "time")
