/* Edit the following line to reflect the full path to your CSV file */
%let csv_file = '/home/u63735368/prostate_cancer/9293SidanaProstateMR_DATA_NOHDRS_2024-03-05_1803.csv';

/* Convert the provider name to uppercase */
data redcap_name;
set redcap;
bx_1_provider = upcase(bx_1_provider);
bx_2_provider = upcase(bx_2_provider);
run;

/* Explore the dataset */

/* Select subjects for whom both bx1_provider and bx2_provider are Dr. Sidana */
data redcap_1;
	set redcap_name;
	where bx_1_provider='SIDANA' and bx_2_provider='SIDANA';
run;

/* Create frequency tables to count total of subjects with and without transition for
Dr. Sidana specifically  */
proc tabulate data=redcap_1;
	class bx1approach bx2approach;
	table bx1approach, bx2approach *N='';
	title 'Number of subjects with and without transition based on Dr. Sidana';
run;

/* Create failure rate table (path) based on three groups (TR to TP, TP
to TP, TR to TR) for Dr. Sidana specifically  */
data temp1;
	set redcap_1;

	if bx1approach = 0 and bx2approach = 1 then
		group=1;

	if bx1approach = 0 and bx2approach = 0 then
		group=2;

	if bx1approach = 1 and bx2approach = 1 then
		group=3;
run;

proc sql;
	create table results1 as select case 
	when group=1 then 'TP-TP' 
	when group=2 then 'TR-TP' 
	when group=3 then 'TR-TR' 
	end as group_label, 
	sum(case when treatment1_pathfail=1 then 1 else 0 end) as fail_count, 
	count(*) as total_count, 
	round(sum(case when treatment1_pathfail=1 then 1 else 0 end)/count(*), 0.01) as fail_rate 
	from temp1 
	group by group;
quit;

proc print data=results1;
	title 'Pathology failure rates based on Dr.Sidana';
run;

/* check if first two bx dates are before outcome dates */
proc sql;
	create table path as select subject_id, bx_1_date, bx_2_date, 
		treatment1_pathfail, treatment1_pathfaildate, bx_1_trg_corespositive, 
		bx_2_trg_corespositive, bx_1_sys_corespositive, bx_2_sys_corespositive, 
		bx1approach, bx2approach, bx_1_provider, bx_2_provider, 
		case when 
		treatment1_pathfaildate ne . and bx_1_date < treatment1_pathfaildate and 
		bx_2_date <=treatment1_pathfaildate then 'met' 
		else 'not met' end as condition 
		from redcap_1;
	select * from path;
	run;

proc sort data=path;
	by condition;
run;

/* create spaghetti plots to see the core positives over time */
data redcap_2;
	set path;

	/* Extract year and month from date columns */
	bx_1_year=year(bx_1_date);
	bx_2_year=year(bx_2_date);
	bx_1_month=month(bx_1_date);
	bx_2_month=month(bx_2_date);
	bx_1_year_month=cats(bx_1_year, put(bx_1_month, z2.));
	bx_2_year_month=cats(bx_2_year, put(bx_2_month, z2.));
	
	array x[2] bx_1_year_month bx_2_year_month;
	array y[2] bx_1_trg_corespositive bx_2_trg_corespositive;
	array z[2] bx_1_sys_corespositive bx_2_sys_corespositive;

	do i=1 to 2;
		x_value=x[i];
		y_value=y[i];
		y2_value=z[i];
		output;
	end;
run;

data redcap_3;
	set redcap_2;
	length group $40;
	where y_value ne .;
	date=input(x_value, yymmn6.);
	format date yymmd7.;

	if treatment1_pathfail=1 and bx1approach ne bx2approach
		then group='fail(transition present)';
	else if treatment1_pathfail=1 and bx1approach = bx2approach 
		then group='fail(no transition)';
	else delete;
run;


proc sgplot data=redcap_3;
	title 'Pathology Failure (targeted): Transition Present vs. No Transition  based on Dr. Sidana';
	styleattrs datacontrastcolors=(bibg brown);
	series x=date y=y_value / group=subject_id markers LINEATTRS=(THICKNESS=2) 
		grouplc=group groupmc=group name='grouping';
	keylegend 'grouping' / type=linecolor;
	xaxis label="year-month";
	yaxis label="# of trg_corespositive" min=-1 max=15;
run;

proc format;
value combined 
0='benign' 1='3+3' 
2='3+4' 3='4+3' 
4='4+4' 5='4+5' 
6='5+4' 7='5+5' 
8='3+5' 9='5+3';
run;

proc format;
value race_fmt 
1 = 'Caucasian'
2 = 'African American'
3, 8, 9 = 'Other'; 
run;

proc format;
value psa_fmt
0 - < 11 = 'psa <= 10'
11 - < 21 = '10 < psa <= 20'
21 - high = 'psa > 20';
run;

proc format;
value volume_fmt
1 - < 41 = 'volume <= 40'
41 - < 61 = '40 < volume <= 60'
61 - high = 'volume > 60';
run;

data redcap_n1;
set redcap_name(keep = subject_id bx1approach bx_1_date bx_1_provider bx_1_sys_gleason bx_1_trg_gleason race ghx_hasfamilyhisotry_v_0 psa_1_value mri_1_prostatevolume);
where bx_1_provider ne '' and bx1approach ne .;
combined = max(bx_1_sys_gleason, bx_1_trg_gleason);
if combined in (0,1) then outcome = 0;
else outcome = 1;
run;

proc print data = redcap_n1(obs = 6);
run;

/* chi-square test */
proc freq data = redcap_n1;
table bx1approach*outcome / chisq;
run;

/* logistic regression */
proc logistic descending data=redcap_n1;
class bx1approach(ref = 'Transrectal')
	  race(ref = 'Caucasian') 
	  ghx_hasfamilyhisotry_v_0(ref = 'No') 
	  psa_1_value(ref = 'psa <= 10') 
	  mri_1_prostatevolume(ref = 'volume <= 40');
	  	  
model outcome = bx1approach race ghx_hasfamilyhisotry_v_0  psa_1_value mri_1_prostatevolume;

format race race_fmt.;
format psa_1_value psa_fmt.;
format mri_1_prostatevolume volume_fmt.;
run;

/* mixed effect logistic regression*/
data redcap_n2;
set redcap_name(keep = subject_id bx1approach bx_1_provider bx_1_sys_gleason bx_1_trg_gleason race ghx_hasfamilyhisotry_v_0 psa_1_value mri_1_prostatevolume);
where bx_1_provider ne '' and bx1approach ne .;
combined = max(bx_1_sys_gleason, bx_1_trg_gleason);
if combined in (0,1) then outcome = 0;
else outcome = 1;
run;

proc glimmix data = redcap_n2;
class bx1approach(ref = 'Transrectal') 
      race(ref = 'Caucasian') 
      ghx_hasfamilyhisotry_v_0(ref = 'No') 
      bx_1_provider
      psa_1_value(ref = 'psa <= 10') 
	  mri_1_prostatevolume(ref = 'volume <= 40');
model outcome(event = '1') = bx1approach race ghx_hasfamilyhisotry_v_0 psa_1_value mri_1_prostatevolume / dist=binary link=logit;
random intercept / subject = bx_1_provider;
format race race_fmt.;
format psa_1_value psa_fmt.;
format mri_1_prostatevolume volume_fmt.;
run;





