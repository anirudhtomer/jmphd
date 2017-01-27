library(nlme)
library(splines)
library(ggplot2)

ticksX = function(max, by){
  scale_x_continuous(breaks = seq(0, max, by = by))
}

ticksY = function(max, by){
  scale_y_continuous(breaks = seq(0, max, by = by))
}

plotRandomProfile = function(count=1, fitted=F){
  pid_sample = sample(x = unique(prias_long$P_ID), size = count)
  plot<-ggplot(data=prias_long[prias_long$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa + 1))) + 
    geom_line(aes(group=P_ID))
  if(fitted==T){
    plot + geom_line(aes(y=fitted, x=visitTimeYears, color=P_ID, group=P_ID)) 
  }
}

ggplot(data=prias_long, aes(x=visitTimeYears, y=log(psa + 1))) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + ticksX(max(prias_long$visitTimeYears), 0.1)

#Since every box plot has 1 entry per person, 
#if the first few box plots have too many time points it means that it is not because
#people came frequently, it could also be due to the fact that people started dropping out later
#and so initially we see more concentration of time points.
#Per person the concentration is not more in the beginning, but overall marginally
#there is more concentration in the beginning
ggplot(data=prias_long[prias_long$P_ID,], aes(factor(visit_number),visitTimeYears)) + 
  geom_boxplot() + stat_summary(fun.data = function(x){
    return(c(y = 1.2, label = length(x))) 
  }, geom = "text", fun.y = median) + ticksY(max(prias_long$visitTimeYears), 0.05)

quantiles_time_years = quantile(prias_long$visitTimeYears, probs = c(0.25, 0.5, 0.75))
# 25%        50%        75% 
# 0.04849315 0.12136986 0.24849315 

#I would take equidistant time points for cutoff
prias_long$logpsa1 = log(prias_long$psa + 1)
prias_long = cbind(prias_long, polyage=poly(prias_long$Age,3))
prias_long = cbind(prias_long, nsFixed=ns(prias_long$visitTimeYears, knots=quantiles_time_years),
                   nsRandom=ns(prias_long$visitTimeYears, knots=quantiles_time_years[1]))
write.csv(prias_long, "prias_long.csv", row.names = F)

psaModel_big = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] +
                 ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85)) * poly(Age,2)[,1] +
                   ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85)) * poly(Age,2)[,2],
                     random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                     data=prias_long[!is.na(prias_long$psa),],
               control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
               method = "ML")

anova(psaModel_big, type="marginal")

psaModel_1 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] +
                     ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85)) * poly(Age,2)[,1],
                   random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                   data=prias_long[!is.na(prias_long$psa),],
                   control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                   method = "ML")

prias_long$fitted[!is.na(prias_long$psa)] = psaModel_1$fitted[,2]
prias_long$residuals[!is.na(prias_long$psa)] = psaModel_1$residuals[,2]

plotRandomProfile(1,T)

psaModel_2 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2)

psaModel_3 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.1, 0.4, 0.7))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2, psaModel_3)

psaModel_4 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30, 0.7))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.03))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2, psaModel_4, psaModel_3 )


psaModel_5 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30, 0.7))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.03, 0.09))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2, psaModel_4, psaModel_3 )




qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
