proc phreg data=prias_id;
model year_discontinued*event_indicator(0) = Age Age*Age/ type3;
run;

/* Cox snell */
proc phreg data=prias_id;
model year_discontinued*event_indicator(0) = Age Age*Age/ type3;
output out = phreg_coxsnell LOGSURV = h;
run;

data phreg_coxsnell;
  set phreg_coxsnell;
  h = -h;  
run;

proc lifetest data=phreg_coxsnell nelson outsurv=survres;
time h*event_indicator(0);
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


/*Longitudinal analysis of PSA*/
proc mixed data=prias_long method=REML;
class P_ID;
model logpsa1 = polyage_1 polyage_2 nsFixed_1-nsFixed_4 
nsFixed_1*polyage_1 nsFixed_2*polyage_1 nsFixed_3*polyage_1 nsFixed_4*polyage_1/solution;
random intercept nsRandom_1 nsRandom_2/type=un subject=P_ID vcorr=1 v=1 g gcorr;
/*random intercept/type=un subject=P_ID vcorr=1 v=1 g gcorr;*/
contrast 'Test main effect splines' nsFixed_1 -1 1, nsFixed_2 -1 1, nsFixed_3 -1 1, nsFixed_4 -1 1;
contrast 'Test interaction effect splines' nsFixed_1*polyage_1 -1 1, nsFixed_2*polyage_1 -1 1, nsFixed_3*polyage_1 -1 1, nsFixed_4*polyage_1 -1 1;
run;
