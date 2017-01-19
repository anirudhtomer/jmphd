library(nlme)
library(splines)
library(ggplot2)

ggplot(data=prias_long, aes(x=visitTimeYears, y=log(psa))) + geom_line(aes(group=P_ID)) + stat_smooth()

pid_sample = sample(x = unique(prias_long$P_ID), size = 500)
pid_sample = unique(prias_long$P_ID)[25]
ggplot(data=prias_long[prias_long$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa))) + geom_line(aes(group=P_ID)) + stat_smooth()


quantiles_time_years = quantile(prias_long$visitTimeYears, probs = c(0.25, 0.5, 0.75))
# 25%        50%        75% 
# 0.04849315 0.12136986 0.24849315 

prias_long$logpsa1 = log(prias_long$psa+1)
prias_long = cbind(prias_long, polyage=poly(prias_long$Age,3))
prias_long = cbind(prias_long, nsFixed=ns(prias_long$visitTimeYears, knots=quantiles_time_years),
                   nsRandom=ns(prias_long$visitTimeYears, knots=quantiles_time_years[1]))
write.csv(prias_long, "prias_long.csv", row.names = F)

psaModel = lme(fixed=log(psa + 1)~Age + 
                 ns(visitTimeYears, knots=c(0.049, 0.121, 0.249)),
                     random = ~ns(visitTimeYears, knots=c(0.049))|P_ID, 
                     data=prias_long[!is.na(prias_long$psa),],
               control = lmeControl(opt = "optim"))
           

anova(psaModel, type="marginal")

qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
