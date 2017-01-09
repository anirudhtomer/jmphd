library(nlme)
library(splines)
library(ggplot2)


ggplot(data=prias_long, aes(x=visitTime, y=log(psa), group=P_ID)) + geom_line()

pid_sample = sample(x = unique(prias_long$P_ID), size = 200)
ggplot(data=prias_long[prias_long$P_ID %in% pid_sample,], aes(x=visitTime, y=log(psa), group=P_ID)) + geom_line()

ggplot(data=prias_long, aes(x=visitTime, y=psa, group=P_ID)) + geom_line() + ylim(0,100)
ggplot(data=prias_long, aes(x=visitTime, y=log(psa))) + geom_point() + stat_smooth()

psaModel = lme(fixed=log(psa)~ns(visitTime, knots=c(33.8, 75.6)) : I(poly(Age,2)[,1]) + 
                 ns(visitTime, knots=c(33.8, 75.6)) : I(poly(Age,2)[,2]) + poly(Age,2)[,1] + poly(Age,2)[,2], 
               random = ~ns(visitTime, knots=c(33.8))|P_ID, 
               data=prias_long[!is.na(prias_long$psa),])

anova(psaModel, type="marginal")

qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
