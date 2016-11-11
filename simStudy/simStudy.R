library(mvtnorm)
library(JMbayes)
library(ggplot2)

orig_b = rmvnorm(n=50, c(0, 0), sigma = matrix(c(40, 10, 10, 4), nrow=2))
time = 1:10
halftime = 1:5

subdata = apply(orig_b, MARGIN = 1, FUN = function(x){
  ret = x[1]
  if(x[1]>0){
    ret = ret + time*x[2] + rnorm(length(time), 0, sd = 0.4)
  }else{
    ret = ret + halftime*x[2] + rnorm(length(halftime), 0, sd = 0.4)                              
  }
  
  ret
})

mydf = data.frame(do.call(rbind, sapply(1:50, FUN = function(x){
  len = length(subdata[[x]])
  dead = if(len==5){1}else{0}
  cbind(rep(x, len),subdata[[x]], 1:len, rep(dead, len))}, 
  simplify = T)))

colnames(mydf) = c("id", "resp", "year","status")

mydf.id = data.frame(do.call(rbind, sapply(1:50, FUN=function(x){
  if(length(subdata[[x]])==5){
    c(x, 5 + runif(n = 1, 1, 5), 1)
  }else{
    c(x, 10 + runif(n=1, 1, 5), 0)
  }
}, simplify = F)))
colnames(mydf.id) = c("id", "year", "status")


ggplot(data = mydf, aes(x=year, y=fittedLme, group=id)) + geom_line()

lmeFit = lme(resp ~ year, data=mydf, random=~year | id)
coxFit = coxph(Surv(mydf.id$year, mydf.id$status) ~ 1, x = TRUE)

jointFit= jointModelBayes(lmeFit, coxFit, 
                                timeVar = "year", n.iter = 30000)

mydf$fittedLme = c(fitted(lmeFit))
mydf$fittedjoint = c(fitted(jointFit, process = "longitudinal", type="subject"))
