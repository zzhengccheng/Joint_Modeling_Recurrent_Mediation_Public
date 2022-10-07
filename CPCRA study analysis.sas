
*Import the dataset from lib ddiddc;
data ddiddc;
set ddiddc.ddiddc;
trt=randgrp-1;
gender=gender-1;
hemobl=hemobl-12;
cd4bl=cd4bl/100;
id=seq;
run;

data one;
set ddiddc;
array event_all{20} MAC PCP1 PCP2 PCP3 CANE1 CANE2 CANE3 CMV WAST 
KSV ADC CRYC TB LYMP PML TOXO CRYS OMYC HIST HZ1;
array time_all{20} T2MAC T2PCP1 T2PCP2 T2PCP3 T2CANE1 T2CANE2 T2CANE3 T2CMV T2WAST 
T2KSV T2ADC T2CRYC T2TB T2LYMP T2PML T2TOXO T2CRYS T2OMYC T2HIST T2HZ1;
aa=1;
do i=1 to 20;
	event=event_all{i};
	stoptime=time_all{i}/30;
	output;
end;
run;
* The dataset for recurrent event;
data two;
set one;
if event=1;
run;

proc sort data=two;
by id stoptime;
run;
* Get the quantiles for recurrent events;
proc univariate data=two noprint;
var stoptime; 
output out=quant_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr; 
run;
data quant_r;
set quant_r;
aa=1;
run;

* The dataset for death event;
data three;
set ddiddc;
by seq;
if first.seq;
stoptime=t2death/30;
event=death*2;		* Set event=2 for death;
aa=1;
run;
* Get the quantiles for death time;
proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
where event=2;
run;
data quant_d;
set quant_d;
aa=1;
run;
* Merge the recurrent and death event times;
data four;
set two three;
run;

proc sort data=four;
by id stoptime;
run;
* Merge data with the quantiles;
data four2;
merge four quant_r quant_d;
by aa;
run;

* Calculate the number of recurrent events as a mediator for death;
data four3;
set four2;
by id;
retain last_stop nevent;
if first.id then do;
	nevent=0;
	start=0;
	stop=stoptime;
	last_stop=stoptime;
end;
else do;
	nevent=nevent+1;
	start=last_stop;
	stop=stoptime;
	last_stop=stoptime;
end;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quant_r {11} qr0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;
array quant_d {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;

array dur_r {10} dur_r1-dur_r10;
array dur_d {10} dur_d1-dur_d10;
array dur_g {10} dur_g1-dur_g10;

array event_r {10} event_r1-event_r10;
array event_d {10} event_d1-event_d10;
array event_g {10} event_g1-event_g10;

array median_r {10} median_r1-median_r10;
array median_d {10} median_d1-median_d10;

do i=1 to 10;
	dur_r{i}=0;
	dur_d{i}=0;
	dur_g{i}=0;
	event_r{i}=0;
	event_d{i}=0;
	
	median_r{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 11;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;
			i=11;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 11;
		if stoptime<=quant_r{i} then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
			i=11;
		end;
		else do;
			dur_r{i-1}=quant_r{i}-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
		end;
	end;

do i=2 to 11;
		if (stoptime<=quant_r{i} and start >quant_r{i-1})then do;
		    dur_g{i-1}=gap;
			end;

        if  (stoptime<=quant_r{i} and start <= quant_r{i-1})then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			if gap >= dur_r{i-1} then dur_g{i-1}=dur_r{i-1};
			else dur_g{i-1}=gap;
			
			i=11;
		end;

		else do;
			dur_g{i-1}=dur_r{i-1};
			
		end;
	end;


	do i=2 to 11;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=11;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;

* For each of mediator "nevent", we calculate the duration in each quantile interval of death time;
data five2;
set five;
array quant {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;
array dur {10} dur1-dur10;

last_start=start;
do i=1 to 10;
	dur{i}=0;
end;

do i=2 to 11;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=11;
	end;
end;
run;	

**************************************************************************************************************
**************************************************************************************************************
************************************************* Model I ****************************************************;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_r1=-3 log_r2=-3 log_r3=-3.2 log_r4=-3 log_r5=-3.2 log_r6=-3 log_r7=-3.3 log_r8=-3.2 log_r9=-3 log_r10=-2.8
	  log_h1=-5.8 log_h2=-5.1 log_h3=-5.3 log_h4=-5 log_h5=-5 log_h6=-4.6 log_h7=-4.6 log_h8=-3.9 log_h9=-3.9 log_h10=-3.5 
	  betax1=0.6 betax2=-.67 betaz=0 betax3=-.1 betax4=-.19 betax5=0
	  etax1=1.46 etax2=-.09 etaz=-.43 etax3=-0.14 etax4=-.54 etax5=0
	  delta1=2 etam=0 log_varc=-1.6;

array dur {10} dur1-dur10;
      * dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array baseh {10} log_h1-log_h10;

base_haz_r=exp(log_r1) * event_r1 + exp(log_r2) * event_r2 + exp(log_r3) * event_r3 
          + exp(log_r4) * event_r4 + exp(log_r5) * event_r5 + exp(log_r6) * event_r6 
          + exp(log_r7) * event_r7 + exp(log_r8)* event_r8 + exp(log_r9) * event_r9 + exp(log_r10) * event_r10;

cum_base_haz_r=exp(log_r1) * dur_r1 + exp(log_r2) * dur_r2 + exp(log_r3) * dur_r3 
              + exp(log_r4) * dur_r4 + exp(log_r5) * dur_r5 + exp(log_r6) * dur_r6 
              + exp(log_r7) * dur_r7 + exp(log_r8)* dur_r8 +exp(log_r9) * dur_r9 + exp(log_r10) * dur_r10;

mu1= betax1 * prevoi + betax2 * gender + betaz*trt + betax3*stratum + betax4*hemobl + betax5* cd4bl + vi; /* for recurrent event */

loglik1=-exp(mu1) * cum_base_haz_r;

sum2=0;
do k=1 to 10;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent);
end;
/* for death event */
mu2= etax1 * prevoi + etax2* gender + etaz*trt + etax3*stratum + etax4 * hemobl + etax5 * cd4bl + delta1 * vi;

loglik0=-exp(mu2) * sum2;

if event=2 then do; /* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3 + exp(log_h4) * event_d4 
              + exp(log_h5) * event_d5 + exp(log_h6) * event_d6 + exp(log_h7) * event_d7 + exp(log_h8)* event_d8 
              + exp(log_h9) * event_d9 + exp(log_h10) * event_d10;
	mu4= etax1 * prevoi + etax2 * gender + etaz * trt + etax3 * stratum + etax4 * hemobl + etax5 * cd4bl 
       + etam * nevent  + delta1 * vi;	
end;

if event=1 then loglik= log(base_haz_r) + mu1 +loglik0; 			/*log likelihood for recurrent event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1 CorrMatParmEst=corr1 CovMatParmEst=cov1; 
run;
*&out is the path for saving results;
data _null_;
	set cov1;
	file "&out.\covI_real.txt";
	put parameter 
		log_r1 log_r2 log_r3 log_r4 log_r5 log_r6 log_r7 log_r8 log_r9 log_r10
		log_h1 log_h2 log_h3 log_h4 log_h5 log_h6 log_h7 log_h8 log_h9 log_h10 
		betax1 betax2 betaz betax3 betax4 betax5
		etax1 etax2 etaz etax3 etax4 etax5
		delta1 etam log_varc;

	run;

data _null_;
	set est1;
	file "&out.\estI_real.txt";
	put Parameter Estimate;
	run;


ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\Estimation setting I.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=est1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Estimation setting I";

	run;
ODS Excel CLOSE;
ODS EXCEL file="&out.\Covariance setting I.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=cov1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Covariance Setting I";

	run;
ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;


**************************************************************************************************************
**************************************************************************************************************
************************************************* Model II ****************************************************;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_r1=-3 log_r2=-3 log_r3=-3.2 log_r4=-3 log_r5=-3.2 log_r6=-3 log_r7=-3.3 log_r8=-3.2 log_r9=-3 log_r10=-2.8
		log_h1=-5.8 log_h2=-5.1 log_h3=-5.3 log_h4=-5 log_h5=-5 log_h6=-4.6 log_h7=-4.6 log_h8=-3.9 log_h9=-3.9 log_h10=-3.5 
		betax1=0.6 betax2=-.67 betaz=0 betax3=-.1 betax4=-.19 betax5=0
		etax1=1.46 etax2=-.09 etaz=-.43 etax3=-0.14 etax4=-.54 etax5=0
		delta1=2 etam=0 delta2=0 log_varc=-1.6;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array dur {10} dur1-dur10;
array baseh {10} log_h1-log_h10;

base_haz_r=exp(log_r1) * event_r1 + exp(log_r2) * event_r2 + exp(log_r3) * event_r3 + exp(log_r4) * event_r4 
          + exp(log_r5) * event_r5 + exp(log_r6) * event_r6 + exp(log_r7) * event_r7 + exp(log_r8)* event_r8 
          + exp(log_r9) * event_r9 + exp(log_r10) * event_r10;

cum_base_haz_r=exp(log_r1) * dur_r1 + exp(log_r2) * dur_r2 + exp(log_r3) * dur_r3 + exp(log_r4) * dur_r4 
              + exp(log_r5) * dur_r5 + exp(log_r6) * dur_r6 + exp(log_r7) * dur_r7 + exp(log_r8) * dur_r8 
              + exp(log_r9) * dur_r9 + exp(log_r10) * dur_r10;

mu1= betax1 * prevoi + betax2* gender + betaz*trt + betax3*stratum + betax4*hemobl + betax5* cd4bl + vi; /* for recurrent event */

loglik1=-exp(mu1) * cum_base_haz_r;

sum2=0;
do k=1 to 10;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + delta2 * nevent * vi);
end;

mu2= etax1 * prevoi + etax2* gender + etaz*trt + etax3*stratum + etax4 * hemobl + etax5 * cd4bl + delta1 * vi;	/* for death event */

loglik0=-exp(mu2) * sum2;

if event=2 then do; /* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3 + exp(log_h4) * event_d4 
              + exp(log_h5) * event_d5 + exp(log_h6) * event_d6 + exp(log_h7) * event_d7 
              + exp(log_h8)* event_d8 + exp(log_h9) * event_d9 + exp(log_h10) * event_d10;
	mu4= etax1 * prevoi + etax2 * gender + etaz * trt + etax3 * stratum + etax4 * hemobl 
       + etax5 * cd4bl + etam * nevent + delta2 * nevent * vi  + delta1 * vi;	
end;

if event=1 then loglik= log(base_haz_r) + mu1 +loglik0; 			/*log likelihood for recurrent event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est2 FitStatistics=fit2 CorrMatParmEst=corr2 CovMatParmEst=cov2; 
run;

data _null_;
	set cov2;
	file "&out.\CovII_real.txt";
	put parameter log_r1 log_r2 log_r3 log_r4 log_r5 log_r6 log_r7 log_r8 log_r9 log_r10
		log_h1 log_h2 log_h3 log_h4 log_h5 log_h6 log_h7 log_h8 log_h9 log_h10 
		betax1 betax2 betaz betax3 betax4 betax5
		etax1 etax2 etaz etax3 etax4 etax5
		delta1 etam delta2 log_varc;

	run;

data _null_;
	set est2;
	file "&out.\EstII_real.txt";
	put Parameter Estimate;
	run;

ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\Estimation setting II.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=est2 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Estimation setting II";

	run;
ODS Excel CLOSE;
ODS EXCEL file="&out.\Covariance setting II.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=cov2 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Covariance Setting II";

	run;

ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;
