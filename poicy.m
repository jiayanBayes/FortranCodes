clear all;

% Define global variables
global REPTITION NOBS XMAT NFC NNC NALT NERROR NFIXED NRANDOM CHOICE ...
       HIGHINC MEDINC DISTANCE EXPRESS OCCUPANCE SOLOEXP HOVEXP2 HOVEXP3 ...
       KE KF HOV2 HOV3 TRANS;  

%%%%%%%%%%% Set Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOMO = 0;

HALTON = 1;

REPTITION = 300;

NOBS = 711;

NALT = 9;

NFC = 16;
NFIXED = 16;

NNC = 6;
NRANDOM = 6;

KE = 4000;
KF = 8000;

%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load hovdata_new.txt;
hovdata = hovdata_new;
CHOICE = hovdata(:,2);
XMAT = hovdata(:,3:12);
DISTANCE = hovdata(:,13);
HIGHINC = hovdata(:,14);
MEDINC = hovdata(:,15);
OCCUPANCE = hovdata(:,16);
clear hovdata;
HIGHINC = HIGHINC + MEDINC;

TRANS   = zeros(NOBS,1);
HOV2    = zeros(NOBS,1);
HOV3    = zeros(NOBS,1);
SOLOEXP = zeros(NOBS,1);
HOVEXP2 = zeros(NOBS,1);
HOVEXP3 = zeros(NOBS,1);
for i = 1:NOBS;
	if CHOICE(i) ~= 1 | CHOICE(i) ~= 4 | CHOICE(i) ~= 7
		TRANS(i) == 1;
	end;
	if CHOICE(i) >= 4 & CHOICE(i) <=6
		HOV2(i) = 1; 
	end;
	if CHOICE(i) >= 7 & CHOICE(i) <= 9
		HOV3(i) = 1; 
	end;
	if CHOICE(i) == 3
		SOLOEXP(i) = 1; 
	end;
	
	if CHOICE(i) == 6
		HOVEXP2(i) = 1; 
	end;
	
	if CHOICE(i) == 9
		HOVEXP3(i) = 1; 
	end;
end;
EXPRESS = SOLOEXP + HOVEXP2 + HOVEXP3;

load beta_revised_het.txt;
b = beta_revised_het;
parameter = [b(1:7);b(9:11);b(13);b(15:16);b(23:25);b(28:31);b(33)];
stop 
%%% New scenario:: cut the HOV share under no-toll by half and
%%% the HOV share is 22.5% under HOV policy%%%
%parameter(11) = parameter(11) + 2.7 * parameter(11);
%parameter(13) = parameter(13) + 1.3 * parameter(13);
%parameter(19) = parameter(19) + 1.8 * parameter(19);

parameter(1:13) = parameter(1:13) * b(34);
parameter(18:20) = parameter(18:20) * b(34);
%%%%%%%%%%% Set the random draws %%%%%%%%%%%%%%%%%%%%%%
PRIME = [2;3;5;7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;79;83;89;97;101];

% Generate Normal Random Numbers Using Halton Sequences
if HALTON == 1
	NERROR = zeros(NOBS, NNC*REPTITION); 
	temp = zeros(NOBS/NALT, NNC*REPTITION);
	for i = 1:NNC;
		gstart = (i-1)*REPTITION + 1;
		gend = (i-1) * REPTITION + REPTITION;
      		seed = PRIME(i);
       		hdraw = halton( (NOBS/NALT)*REPTITION, 0, seed, 1);
       		hdraw = transpose(hdraw);
       		for j = 1:(NOBS/NALT);
       			temp(j,gstart:gend) =...
            		hdraw( 1, ((j-1)*REPTITION + 1) : ((j-1)*REPTITION+REPTITION) );
       		end;
		clear hdraw;
	end;
	clear gstart gend;
	
	ini=0;
      	for i = 1:(NOBS/NALT);
       		gstart = ini+1;
       		gend = ini+NALT;
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(NALT,1) );
       		ini=ini+NALT;
       	end;
       	clear temp;
   
else
	NERROR = zeros(NOBS, NNC*REPTITION);
    	temp = randn((NOBS/NALT), NNC*REPTITION);
    	ini=0;
      	for i = 1:(NOBS/NALT);
       		gstart = ini+1;
       		gend = ini+NALT;
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(NALT,1) );
       		ini=ini+NALT;
       	end;
       	clear temp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Congested case with het. %%%
ntotal = 17570;
cons = 12.65;

rho = 0.2;
%%%%%% Homogeneous Case %%%%%%%%%
if HOMO == 1
	tv1 = mean(DISTANCE);
	tv2 = mean(HIGHINC);
	
	for i = 1:length(DISTANCE);
		DISTANCE(i,1) = tv1;
		HIGHINC(i,1) = tv2;
	end;
	
	parameter(17) = 0;
	parameter(21) = 0;
	
	%% We can also eliminate the systematic preferences for express lane and transponder%%%
	parameter(7) = parameter(7) + 0.2532 * parameter(8) + parameter(9); 
	parameter(8:9) = 0;
	
	parameter(10) = parameter(10) + 0.3671 * parameter(14) + 0.6203 * parameter(15) + 2.9747 * parameter(16);
	parameter(14:16) = 0;
end;

%%%% equilibrium in base case: notoll, and nocarpool%%%%%% 

[junk, pre] = equilibrium(ntotal, cons, rho, 3.30, 0, parameter, 0);
veupdate = sum(junk(:,3) + junk(:,6)/2 + junk(:,9)/3);
vfupdate = sum(junk(:,1) + junk(:,2) + junk(:,4)./2 + junk(:,5)./2 + junk(:,7)./3 + junk(:,8)./3);

neupdate = sum(junk(:,3) + junk(:,6)/1 + junk(:,9)/1);
nfupdate = sum(junk(:,1) + junk(:,2) + junk(:,4)./1 + junk(:,5)./1 + junk(:,7)./1 + junk(:,8)./1);

pfupdate = sum(junk(:,4) + junk(:,5) + junk(:,6)./1 + junk(:,7)./1 + junk(:,8)./1 + junk(:,9)./1);
shov = pfupdate / ntotal
pfupdate = sum(junk(:,7)./1 + junk(:,8)./1 + junk(:,9)./1);
shov3 = pfupdate / ntotal

tenotoll = supply(KE, veupdate)
tfnotoll = supply(KF, vfupdate)
renotoll = (tenotoll - 9.23) * 0.3785;
rfnotoll = (tfnotoll - 9.23) * 0.3785;

[pa, pn1] = demand (ntotal,cons, rho, parameter, 3.30*1.1, 0*1.1, tenotoll*1.1, tfnotoll*1.1, renotoll*1.1, rfnotoll*1.1, 0);

vnew = sum(pn1(:,3) + pn1(:,6)/2 + pn1(:,9)/3);

(vnew - veupdate) / veupdate;

elas = ( sum(sum(pn1)) - sum(sum(junk)) ) / sum(sum(junk))

STOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tini = 0;

dplot = zeros(100, 3);
for i = 1: 100;
	i
	delt = (i - 1) * 0.05;
	tao = tini + delt;
	
	tenotoll = 20;
	tfnotoll = 20;
	renotoll = (tenotoll - 9.23) * 0.3785;
	rfnotoll = (tfnotoll - 9.23) * 0.3785;

	dplot(i,1) = tao;

	[pa, pn1] = demand (ntotal1,cons1, rho, parameter, tao, tao, tenotoll, tfnotoll, renotoll, rfnotoll, 6);
	[pa, pn2] = demand (ntotal2,cons2, rho, parameter, tao, tao, tenotoll, tfnotoll, renotoll, rfnotoll, 6);
	
	veupdate = sum(pn1(:,3) + pn1(:,6)/2 + pn1(:,9)/3);
	vfupdate = sum(pn1(:,1) + pn1(:,2) + pn1(:,4)./2 + pn1(:,5)./2 + pn1(:,7)./3 + pn1(:,8)./3);
	dplot(i,2) = veupdate + vfupdate;

	veupdate = sum(pn2(:,3) + pn2(:,6)/2 + pn2(:,9)/3);
	vfupdate = sum(pn2(:,1) + pn2(:,2) + pn2(:,4)./2 + pn2(:,5)./2 + pn2(:,7)./3 + pn2(:,8)./3);
	dplot(i,3) = veupdate + vfupdate;
end
		
