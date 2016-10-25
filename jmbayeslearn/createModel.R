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

pbc2$status2 = as.numeric(pbc2$status!="alive")
pbc2.id$status2 = as.numeric(pbc2.id$status != "alive")

lmeFit.pbc1 = lme(log(serBilir) ~ ns(year,2), data=pbc2, random=~ns(year,2) | id)
coxFit.pbc1 = coxph(Surv(years, status2) ~ drug * age, data = pbc2.id, x = TRUE)

jointFit.pbc1 = jointModelBayes(lmeFit.pbc1, coxFit.pbc1, 
                                timeVar = "year", n.iter = 30000)

dForm = list(fixed = ~0 + dns(year, 2), random = ~0 + dns(year, 2), 
              indFixed = 2:3, indRandom = 2:3)

jointFit.pbc12 = update(jointFit.pbc1, param="td-both", extraForm=dForm)

tf1 = function(x, data) {
  cbind(x, "^2" = x * x)
}

tf2 = function(x, data) {
  cbind(x, "D-penicil" = x * (data$drug == "D-penicil"))
}

jointFit.pbc15 = update(jointFit.pbc12, transFun = list("value" = tf1, "extra" = tf2))


