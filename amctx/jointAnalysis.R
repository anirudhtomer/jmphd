################
# Patient 346 is not having any PCR measurements...so exclude him
################



library(JMbayes)

jointFit_ns1random = jointModelBayes(longModel_ns2_slope_yr_mv, coxModel, timeVar = "tx_s_years", n.iter = 30000)
dForm <- list(fixed = ~ 0 + dns(tx_s_years, knots=c(200, 500, 1000)/365, Boundary.knots = c(0, 11)), 
              random = ~ 0 + dns(tx_s_years, knots = c(200, 500)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4)

jointFit_asso_slope <- update(jointFit_ns1random, param = "td-both", extraForm = dForm)

qplot(y=residuals(longModel, type = "response"), x=fitted(longModel), geom=c("point", "smooth"))
qplot(y=residuals(jointFit, process="Longitudinal", type = "marginal"), x=fitted(jointFit, process="Longitudinal", type = "marginal"), geom=c("point", "smooth"))

plot(density(jointFit_ns1random$mcmc$betas[,2]))


#########################################
# Using mv joint model function
###########################################
mvJointFit_asso_tdvalue=mvJointModelBayes(longModel_ns2_slope_yr_mv_both, coxModel, timeVar = "tx_s_years")

Forms <- list("value" = "value",
              "value" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(200, 500, 1000)/365, Boundary.knots = c(0, 11)), 
                             random = ~ 0 + dns(tx_s_years, knots = c(200, 500)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                             name = "slope"))

mvJointFit_asso_tdboth <- update(mvJointFit_asso_tdvalue, Formulas = Forms)
