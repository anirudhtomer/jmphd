library(lme4)
library(splines)
library(ggplot2)

quantiles_time_years = round(quantile(prias_long[!is.na(prias_long$gleason),]$visitTimeYears, probs = seq(0,1, by=0.005)), 3)
quantiles_time_years = unique(quantiles_time_years)

odds = vector("numeric", length(quantiles_time_years))
quantiles_time_years = c(-1, quantiles_time_years)
for(i in 2:length(quantiles_time_years)){
  prias_df_temp = prias_long[prias_long$visitTimeYears<=quantiles_time_years[i] & 
                               prias_long$visitTimeYears > quantiles_time_years[i-1],]
  odds[i-1] = table(prias_df_temp$digleason)["Low"]/table(prias_df_temp$digleason)["High"]
}

qplot(y=log(odds), x=quantiles_time_years[-1], geom=c("point","line", "smooth")) + ticksX(1.5, 0.05)

ggplot(data = prias_long, aes(y=as.numeric(digleason),x=visitTimeYears)) + geom_point() + stat_smooth()

#Model fitting
glmer(digleason ~ visitTimeDays + (1|P_ID),
        data=prias_long[!is.na(prias_long$digleason),], family="binomial")


model_gleasons = mvglmer(list(digleason ~Age + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) + 
                       (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05))|P_ID)),
        data = temp, families = list(binomial))

model_dre = mvglmer(list(didre ~ Age + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) + 
                       (ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5))|P_ID)),
                data = prias_long, families = list(binomial))

#subject 957 doesn't have any gleason. removing it

model_all = mvglmer(list(

                    didre ~ Age + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) +
                       (ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5))|P_ID),

                    digleason ~ Age + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) +
                           (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05))|P_ID),
                    
                    logpsa1 ~ Age  +
                      ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85), Boundary.knots = c(0,1.5)) * Age+
                      (ns(visitTimeYears, knots=c(0.1), Boundary.knots = c(0,1.5))|P_ID)

                    ),data = prias_long[prias_long$P_ID!=957,], 
                    families = list(binomial, binomial, gaussian))



mvglmer(list(digleason ~ Age + visitTimeYears +
               d_age +  tx_dgf + is_cni + d_bmi + tx_hla + tx_previoustx + 
               tx_pra + tx_cit + tx_dial_days + tx_dm + rec_bmi + 
               ns(tx_s_years, knots=(c(100, 300, 1000)/365), Boundary.knots = c(0, 11)) + 
               (ns(tx_s_years, knots=(c(100, 300)/365), Boundary.knots = c(0, 11))|amctx)),
        data = prias_long[!is.na(prias_long$digleason),], families = list(binomial))


glmmPQL(digleason ~  visitTimeYears, random = ~ visitTimeYears | P_ID,
        family = "binomial", data = prias_long[!is.na(prias_long$digleason),])


anova(model_dre_1, type="marginal")

qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
