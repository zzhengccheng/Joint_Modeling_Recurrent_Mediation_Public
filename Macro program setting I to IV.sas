
***********************************************************************************
*************************Setting I***************************************************;

%macro jointmodel1(Data=,mod=);

%do ii=1 %to 200;
title "Iteration &ii.";
*Import the dataset;
*&InD is the path for import the datasets;
PROC IMPORT OUT= WORK.sim1 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
*Import the jump points;
data quant_r;
infile "&InD.\quant.txt";
input qr0 qr2 qr4 qr6 qr10 aa;
run;

data quant_d;
infile "&InD.\quant.txt";
input qd0 qd2 qd4 qd6 qd10 aa;
run;

data four;
set sim1;
aa=1;
run;

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
array quant_r {5} qr0 qr2 qr4 qr6 qr10;
array quant_d {5} qd0 qd2 qd4 qd6 qd10;

array dur_r {4} dur_r1-dur_r4;
array dur_d {4} dur_d1-dur_d4;

array event_r {4} event_r1-event_r4;
array event_d {4} event_d1-event_d4;

array median_r {4} median_r1-median_r4;
array median_d {4} median_d1-median_d4;

do i=1 to 4;
	dur_r{i}=0;
	dur_d{i}=0;
	event_r{i}=0;
	event_d{i}=0;
	median_r{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 5;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;
			i=5;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 5;
		if stoptime<=quant_r{i} then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
			i=5;
		end;
		else do;
			dur_r{i-1}=quant_r{i}-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
		end;
	end;

	do i=2 to 5;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=5;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;

* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;
data five2;
set five;
array quant {5} qd0 qd2 qd4 qd6 qd10;
array dur {4} dur1-dur4;

last_start=start;
do i=1 to 4;
	dur{i}=0;
end;

do i=2 to 5;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=5;
	end;
end;
run;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_lam_m1=-1.2 log_lam_m2=-0.7 log_lam_m3=-0.5 log_lam_m4=-0.43
		log_lam_y1=-3.4 log_lam_y2=-2.5 log_lam_y3=-2.0 log_lam_y4=-0.2
		betax=0.2 betaz=0.35
		etaz=0.35 etax=0.15
		etam=0.25
        delta1=1
		log_varc=-0.7;

* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array dur {4} dur1-dur4;
array baseh {4} log_lam_y1-log_lam_y4;

base_haz_r=exp(log_lam_m1) * event_r1 + exp(log_lam_m2)* event_r2 + exp(log_lam_m3) * event_r3 + exp(log_lam_m4) * event_r4;

cum_base_haz_r=exp(log_lam_m1) * dur_r1 + exp(log_lam_m2) * dur_r2 + exp(log_lam_m3) * dur_r3 + exp(log_lam_m4) * dur_r4;

mu1= betaz * X1 + betax * X2 + vi;	/* for recurrent event */

loglik1=-exp(mu1) * cum_base_haz_r;

sum2=0;
do k=1 to 4 ;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent);
	
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;	/* for death event */

loglik0=-exp(mu2) * sum2;


if event=2 then do;  /* for death event */
    base_haz_d=exp(log_lam_y1) * event_d1 + exp(log_lam_y2) * event_d2 + exp(log_lam_y3) * event_d3 + exp(log_lam_y4) * event_d4;
	
	mu4= etaz * X1 + etax * X2 + etam * nevent  + delta1 * vi;	
   
end;

if event=1 then loglik= log(base_haz_r) + mu1 +loglik0; 			/*log likelihood for recurrent event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1;							/*log likelihood for censoring */

model id ~ general(loglik);

random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est&mod&ii FitStatistics=fit&mod&ii CorrMatParmEst=corr&mod&ii CovMatParmEst=cov&mod&ii; 
run;

%end;

%mend;

%jointmodel1(Data=sim1_alldata,mod=1);

/* &out is the saving file path*/
*Output the covariance;

%macro outcov(mod=,filename=,delta2=);
%do ii=1 %to 200;
data _null_;
	set cov&mod&ii;
	file "&out.\&filename&ii..txt";
	put parameter log_lam_m1 log_lam_m2 log_lam_m3 log_lam_m4
		log_lam_y1 log_lam_y2 log_lam_y3 log_lam_y4
		betax betaz
		etaz etax
        etam
        delta1 &delta2
		log_varc;
	run;
%end;
%mend;

*Output the estimation;
%macro outest(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set est&mod&ii;
	file "&out.\&filename&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outcov(mod=1,filename=covIlog,delta2=);
%outest(mod=1,filename=estIlog);





***********************************************************************************
*************************Setting II***************************************************;


%macro jointmodel2(Data=,mod=);

%do ii=1 %to 200;
title "Iteration &ii.";
*&InD is the path for import the datasets;
PROC IMPORT OUT= WORK.sim2 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data quant_r;
infile "&InD.\quant.txt";
input qr0 qr2 qr4 qr6 qr10 aa;
run;

data quant_d;
infile "&InD.\quant.txt";
input qd0 qd2 qd4 qd6 qd10 aa;
run;

data four;
set sim2;
aa=1;
run;

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
array quant_r {5} qr0 qr2 qr4 qr6 qr10;
array quant_d {5} qd0 qd2 qd4 qd6 qd10;

array dur_r {4} dur_r1-dur_r4;
array dur_d {4} dur_d1-dur_d4;

array event_r {4} event_r1-event_r4;
array event_d {4} event_d1-event_d4;

array median_r {4} median_r1-median_r4;
array median_d {4} median_d1-median_d4;

do i=1 to 4;
	dur_r{i}=0;
	dur_d{i}=0;
	event_r{i}=0;
	event_d{i}=0;
	median_r{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 5;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;
			i=5;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 5;
		if stoptime<=quant_r{i} then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
			i=5;
		end;
		else do;
			dur_r{i-1}=quant_r{i}-quant_r{i-1};
			median_r{i-1}=quant_r{i-1}+dur_r{i-1}/2; /* Get the median of each interval */
			lastmed_r=median_r{i-1};
		end;
	end;

	do i=2 to 5;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=5;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;

* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;
data five2;
set five;
array quant {5} qd0 qd2 qd4 qd6 qd10;
array dur {4} dur1-dur4;

last_start=start;
do i=1 to 4;
	dur{i}=0;
end;

do i=2 to 5;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=5;
	end;
end;
run;	

* Model II in Biometrics paper;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log_lam_m1=-1.2 log_lam_m2=-0.7 log_lam_m3=-0.5 log_lam_m4=-0.43
		log_lam_y1=-3.4 log_lam_y2=-2.5 log_lam_y3=-2.0 log_lam_y4=-0.2
		betax=0.2 betaz=0.35
		etaz=0.35 etax=0.15
		etam=0.25
        delta1=1 delta2=-0.5
		log_varc=-0.7;

* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
array dur {4} dur1-dur4;
array baseh {4} log_lam_y1-log_lam_y4;


base_haz_r=exp(log_lam_m1) * event_r1 + exp(log_lam_m2)* event_r2 + exp(log_lam_m3) * event_r3 + exp(log_lam_m4) * event_r4;

cum_base_haz_r=exp(log_lam_m1) * dur_r1 + exp(log_lam_m2) * dur_r2 + exp(log_lam_m3) * dur_r3 + exp(log_lam_m4) * dur_r4;

mu1= betaz * X1 + betax * X2 + vi;			/* for recurrent event */

loglik1=-exp(mu1) * cum_base_haz_r;


sum2=0;
do k=1 to 4;
	/* cumulative baseline hazard for time dependent measure */
sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + delta2 * nevent * vi);
	
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;	/* for death event */
loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_lam_y1) * event_d1 + exp(log_lam_y2) * event_d2 + exp(log_lam_y3) * event_d3 + exp(log_lam_y4) * event_d4;
	mu4= etaz * X1 + etax * X2 + etam * nevent + delta1 * vi + delta2 * nevent * vi;
	
end;

if event=1 then loglik=log(base_haz_r) + mu1 +loglik0; 			/*log likelihood for recurrent event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est&mod&ii FitStatistics=fit&mod&ii CorrMatParmEst=corr&mod&ii CovMatParmEst=cov&mod&ii;

run;
%end;

%mend;


%jointmodel2(Data=sim2_alldata,mod=2);
*Export the estimation and covariance results;
%outcov(mod=2,filename=covIIlog,delta2=delta2);
%outest(mod=2,filename=estIIlog);



***********************************************************************************
*************************Setting III***************************************************;

*Datasets sim3_alldata generating using log gamma distribution;
%jointmodel2(Data=sim3_alldata,mod=3);
*Export the estimation and covariance results;
%outcov(mod=3,filename=covIIIlog,delta2=delta2);
%outest(mod=3,filename=estIIIlog);


**************************************************************************************
*************************Setting IV***************************************************;

*Fit model without the interaction term while the data are generated with the interaction term;
%jointmodel1(Data=sim2_alldata,mod=4);
*Export the estimation and covariance results;
%outcov(mod=4,filename=covIVlog,delta2=);
%outest(mod=4,filename=estIVlog);
