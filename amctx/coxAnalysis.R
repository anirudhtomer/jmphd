library(survival)
library(penalized)
library(glmnet)
library(doParallel)

colsOfInterest = c("rec_age", "rec_gender",
                   "tx_previoustx", "d_age", "d_gender", "d_bmi", "rec_bmi",
                   "tx_hla", "tx_pra", "tx_dgf", "tx_cit", "is_nr", "is_pred", 
                   "is_aza", "is_cni", "is_mmf", "is_mtor", "ah_nr", "ah_diur",      
                   "ah_ace", "ah_arb", "ah_raasi", "ah_bb", "ah_ccb", "statin", 
                   "dm_oad", "dm_insulin", "tx_dm","tx_hvdis", "rr_sys",
                   "rr_dias", "rr_map", "tx_dial_days", "d_type", "d_cadaveric")

##################
#Step 1: simple additive model with all covariates
##################
opt <- optL1(Surv(days_tx_gl, gl_failure), 
             penalized =~rec_age + rec_gender +
               tx_previoustx + d_age + d_gender + d_bmi +
               rec_bmi + tx_hla + tx_pra + tx_dgf + 
               tx_cit+ is_nr + is_aza + 
               is_cni + is_mmf  + ah_nr + 
               ah_diur + ah_ace + ah_arb + 
               ah_raasi + ah_bb + ah_ccb  + 
               dm_oad + dm_insulin + tx_dm + tx_hvdis + 
               rr_sys + rr_dias + rr_map + 
               tx_dial_days + d_type + d_cadaveric, 
             data = amctx.id, fold = 10, standardize = T, minlambda1 = 30, maxiter = 100)

fit1 <- profL1(Surv(days_tx_gl, gl_failure), 
               penalized =~ rec_age + rec_gender +
                 tx_previoustx + d_age + d_gender + d_bmi +
                 rec_bmi + tx_hla + tx_pra + tx_dgf + 
                 tx_cit+ is_nr + is_aza + 
                 is_cni + is_mmf  + ah_nr + 
                 ah_diur + ah_ace + ah_arb + 
                 ah_raasi + ah_bb + ah_ccb  + 
                 dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                 rr_sys + rr_dias + rr_map + 
                 tx_dial_days + d_type + d_cadaveric, 
               data = amctx.id, standardize = T,fold=10, plot=TRUE, minlambda1 = 1, maxlambda1 = 30)

penalized_model = penalized(Surv(days_tx_gl, gl_failure), 
                            penalized = ~rec_age + rec_gender +
                              tx_previoustx + d_age + d_gender + d_bmi +
                              rec_bmi + tx_hla + tx_pra + tx_dgf + 
                              tx_cit+ is_nr + is_aza + 
                              is_cni + is_mmf  + ah_nr + 
                              ah_diur + ah_ace + ah_arb + 
                              ah_raasi + ah_bb + ah_ccb  + 
                              dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                              rr_sys + rr_dias + rr_map + 
                              tx_dial_days + d_type + d_cadaveric,
                            data = amctx.id, lambda1=10, standardize=TRUE)

coeff = slot(penalized_model, "penalized")
attributes(coeff)$names[coeff != 0]

model1 = coxph(Surv(days_tx_gl, gl_failure) ~ 1+rec_bmi + d_age + d_bmi + d_type + tx_cit + 
        tx_previoustx + tx_hla + rec_age + tx_dial_days + tx_dm + tx_pra,
      data = amctx.id)

zeroVar = c("statin", "is_pred", "is_mtor")
#What to do with ah_nr(1 col with 0), dm_insulin(infinite beta for no reason), d_cadaveric(singular)

modelNull = coxph(Surv(days_tx_gl, gl_failure) ~ 1,
                  data = amctx.id)

model1 = coxph(Surv(days_tx_gl, gl_failure) ~ rec_age  + rec_gender +
                 tx_previoustx + d_age + d_gender + d_bmi +
                 rec_bmi + tx_hla + tx_pra + tx_dgf + 
                 tx_cit+ is_nr  + is_aza + 
                 is_cni + is_mmf   + 
                 ah_diur + ah_ace + ah_arb + 
                 ah_raasi + ah_bb + ah_ccb  + 
                dm_insulin + tx_dm + tx_hvdis + 
                 rr_sys + rr_dias + rr_map + 
                 tx_dial_days + d_type,
               data = amctx.id)
colnames(model.matrix(model1))

#All 3 lead to the same conclusion
stepAIC(model1, direction = "backward")
stepAIC(modelNull,direction="forward",scope=list(upper=model1,lower=modelNull))
stepAIC(modelNull,direction="both",scope=list(upper=model1,lower=modelNull))

########################################################
#AIC + results from shrinkage + what clinicians recommend
#######################################################
opt <- optL1(Surv(days_tx_gl, gl_failure), 
             penalized =~d_age + d_bmi + d_type + tx_cit +
               tx_hla + rec_age + tx_previoustx + tx_dial_days + tx_dm + rec_bmi +
               tx_pra + is_mmf + ah_diur + d_gender, 
             data = amctx.id, fold = 10, standardize = T, minlambda1 = 30, maxiter = 100)

fit1 <- profL1(Surv(days_tx_gl, gl_failure), 
               penalized =~ d_age + d_bmi + d_type + tx_cit +
                 tx_hla + rec_age + tx_previoustx + tx_dial_days + tx_dm + rec_bmi +
                 tx_pra + is_mmf + ah_diur + d_gender, 
               data = amctx.id, standardize = T,fold=10, plot=TRUE, minlambda1 = 1, maxlambda1 = 30)

penalized_model = penalized(Surv(days_tx_gl, gl_failure), 
                            penalized = ~d_age + d_bmi + d_type + tx_cit +
                              tx_hla + rec_age + tx_previoustx + tx_dial_days + tx_dm + rec_bmi +
                              tx_pra + is_mmf + ah_diur + d_gender,
                            data = amctx.id, lambda1=10, standardize=TRUE)
coeff = slot(penalized_model, "penalized")
attributes(coeff)$names[coeff != 0]

coxModel = coxph(Surv(years_tx_gl, gl_failure) ~ d_age + d_bmi + d_type,
               data = amctx.id, x=T, model=T)


plot(survfit(Surv(years_tx_gl, gl_failure)~1, data=amctx.id))


time.dep.zph <- cox.zph(coxModel)
for(i in 1:17){
  plot(time.dep.zph[i])
  abline(h = 0, lty=3)
}

coxSnellRes=-log(amctx.id$gl_failure-resid(coxModel,type="martingale"))
fitres=survfit(coxph(Surv(coxSnellRes,amctx.id$gl_failure)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function')
abline(0,1,col='red',lty=2)
