/*coxModel = coxph(Surv(days_tx_gl, gl_failure) ~ d_age + d_bmi + d_type + tx_cit +
                              tx_hla + rec_age + tx_previoustx + tx_dial_days + tx_dm + rec_bmi +
                              tx_pra + is_mmf + ah_diur + d_gender,
               data = amctx.id, x=T)*/
data amctx;
set amctx;
logdays_tx_gl = log(days_tx_gl);
run;

/* Cox snell */
proc phreg data=amctx;
class d_type tx_hla tx_previoustx tx_dm;
model days_tx_gl*gl_failure(0) = d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx tx_dial_days tx_dm rec_bmi tx_pra;
output out = phreg_coxsnell LOGSURV = h;
run;

proc phreg data=amctx;
class d_type tx_hla tx_previoustx tx_dm;
model days_tx_gl*gl_failure(0) = d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx tx_dial_days tx_dm rec_bmi tx_pra
d_age*logdays_tx_gl d_bmi*logdays_tx_gl rec_age*days_tx_gl;
output out = phreg_coxsnell LOGSURV = h;
run;

data phreg_coxsnell;
  set phreg_coxsnell;
  h = -h;  
run;

proc lifetest data=phreg_coxsnell nelson outsurv=survres;
time h*gl_failure(0);
run;

data survres;
set survres;
ls=(-log(survival));
run;

proc sgplot data=survres NOAUTOLEGEND;
step x=h y=ls;
xaxis label="Cox Snell residuals";
yaxis label="Estimated cumulative hazard rate";
lineparm x=0 y=0 slope=1 /lineattrs=(color=black);
run;

/* Schoenfeld residuals */
proc phreg data=amctx;
class d_type tx_hla tx_previoustx tx_dm;
model days_tx_gl*gl_failure(0) = d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx tx_dial_days tx_dm rec_bmi tx_pra;
output out = phreg_schoenfeld wtressch = d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx tx_dial_days tx_dm rec_bmi tx_pra;
run;

proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = d_age / clm;
  yaxis label="Standardize Schoenfeld residuals: d_Age";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;   
  loess x = days_tx_gl y = d_bmi / clm;
  yaxis label="Standardized Schoenfeld residuals: d_bmi";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = tx_cit / clm;
  yaxis label="Standardized Schoenfeld residuals: tx_cit";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = rec_age / clm;
  yaxis label="Standardized Schoenfeld residuals: rec_age";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = tx_dial_days / clm;
  yaxis label="Standardized Schoenfeld residuals: tx_dial_days";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = rec_bmi / clm;
  yaxis label="Standardized Schoenfeld residuals: rec_bmi";
run;
proc sgplot data=phreg_schoenfeld NOAUTOLEGEND;  
  loess x = days_tx_gl y = tx_pra / clm;
  yaxis label="Standardized Schoenfeld residuals: tx_pra";
run;


/* Martingale residuals age */
proc phreg data=amctx;
class d_type tx_hla tx_previoustx tx_dm;
model days_tx_gl*gl_failure(0) = d_age d_bmi d_type tx_cit tx_hla rec_age tx_previoustx  tx_dm rec_bmi tx_pra;
output out = phreg_martingale resmart=mres;
run;

proc sgplot data=phreg_martingale NOAUTOLEGEND;  
  loess x = tx_dial_days y = mres / clm;
  yaxis label="Martingale Residuals: Age left out";
run;

/* Martingale residuals prio */
proc phreg data=convict;
class fin mar wexp race paro educ;
model week*arrest(0) = fin mar wexp age;
output out = phreg_martingale resmart=mres;
run;

proc sgplot data=phreg_martingale NOAUTOLEGEND;  
  loess x = prio y = mres / clm;
  yaxis label="Martingale Residuals: Prior convictions left out";
run;

/* DFBETA and LD */
proc phreg data=convict_nolabel;
class fin mar wexp race paro educ;
model week*arrest(0) = fin mar mcage prio wexp mar*wexp mcage*fin;
output out=dfbetald dfbeta=fin mar mcage prio wexp mar_wexp mcage_fin ld=ld;
run;

data dfbetald;
set dfbetald;
seqno=_n_;
if(ld>0.10) then obsld=_n_;else obsld=.;
if(mar_wexp>= 0.1 or mar_wexp<=-0.1) then obsmarwexpno=_n_; else obsmarwexpno=.;
if(mcage_fin>= 0.01 or mcage_fin<=-0.01) then obsmcagefinno=_n_; else obsmcagefinno=.;
if(prio2>= 0.005 or prio2<=-0.005) then obspriono=_n_; else obspriono=.;
if(mar2<=-0.2) then obsmarno=_n_; else obsmarno=.;
if(fin2>= 0.025 or fin2<=-0.025) then obsfinno=_n_; else obsfinno=.;
if(mcage2>= 0.01) then obsageno=_n_; else obsageno=.;
if(wexp2>= 0.1 or wexp2<=-0.1) then obswexpno=_n_; else obswexpno=.;
run;

proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=ld/datalabel=obsld;
xaxis label="Observaton Number";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=mar_wexp/datalabel=obsmarwexpno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Marriage * Work Experience";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=mcage_fin/datalabel=obsmcagefinno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Age * Financial Aid";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=prio2/datalabel=obspriono;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Prior convictions";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=fin2/datalabel=obsfinno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Financial Aid";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=mar2/datalabel=obsmarno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Marriage Status";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=mcage2/datalabel=obsageno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Age";
run;
proc sgplot data=dfbetald NOAUTOLEGEND;
scatter x=seqno y=wexp2/datalabel=obswexpno;
xaxis label="Observaton Number";
yaxis label="Difference in parameter for Work experience";
run;

/* without 427,292, 119*/
data convict_noinfluence;
set convict_nolabel;
run;
proc sql;
delete from convict_noinfluence where _=427 or _=292 or _=119;
run;

proc phreg data=convict_noinfluence;
class fin mar wexp race paro educ;
model week*arrest(0) = fin mar mcage prio wexp mar*wexp mcage*fin;
run;


/*model selection*/
proc phreg data=convict_nolabel;
class fin mar wexp race paro educ;
model week*arrest(0) = fin mar mcage educ mcprio wexp race paro 
employed agetime wexptime
fin * mar
fin * mcage
fin * educ
fin * mcprio
fin * wexp
fin * race
fin * paro
mar * mar
mar * mcage
mar * educ
mar * mcprio
mar * wexp
mar * race
mar * paro
mcage * mcage
mcage * educ
mcage * mcprio
mcage * wexp
mcage * race
mcage * paro
educ * educ
educ * mcprio
educ * wexp
educ * race
educ * paro
mcprio * mcprio
mcprio * wexp
mcprio * race
mcprio * paro
wexp * wexp
wexp * race
wexp * paro
race * race
race * paro
paro * paro
/selection=stepwise hier=multiple;
array emp{*} emp1-emp52;
do i=1 to 52;
if(week>=4) then employed = (emp[week]+emp[week-1]+emp[week-2]+emp[week-3])/4; /* or emp[week]*/
else if(week=3) then employed = (emp[week]+emp[week-1]+emp[week-2])/3;
else if(week=2) then employed = (emp[week]+emp[week-1])/2;
else if(week=1) then employed = emp[week];
end;
employed=employed*100;
agetime=mcage*week;
wexptime=wexp*week;
run;

/* Final model */
proc phreg data=convict_nolabel;
class fin mar wexp race paro educ;
model week*arrest(0) = fin mar mcage mcprio wexp 
employed agetime wexptime mar*wexp fin*mcage;
array emp{*} emp1-emp52;
do i=1 to 52;
if(week>=4) then employed = (emp[week]+emp[week-1]+emp[week-2]+emp[week-3])/4; /* or emp[week]*/
else if(week=3) then employed = (emp[week]+emp[week-1]+emp[week-2])/3;
else if(week=2) then employed = (emp[week]+emp[week-1])/2;
else if(week=1) then employed = emp[week];
end;
employed=employed*100;
agetime=mcage*week;
wexptime=wexp*week;
run;

/*4a, 4b*/
proc lifereg data=convict;
class fin wexp;
model week*arrest(0) = fin wexp fin*wexp/distribution=weibull covb;
run;

proc phreg data=convict;
class fin mar wexp race paro educ;
model week*arrest(0) =fin fintime;
fintime=fin*log(week);
run;
