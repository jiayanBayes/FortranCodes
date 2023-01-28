clear all;

% Define global variables
global REPTITION NOBS XMAT NFC NNC NALT NERROR NFIXED NRANDOM CHOICE ...
       HIGHINC MEDINC DISTANCE IDIOSYNCRATIC EXPRESS OCCUPANCE SOLOEXP HOVEXP2 HOVEXP3...
       NPOP CONS PARAMETER KE KF HOV2 HOV3 TRANS;  

%%%%%%%%%%% Set Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOMO = 0;

HALTON = 1;

REPTITION = 100;

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
parameter(11) = parameter(11) + 2.7 * parameter(11);
parameter(13) = parameter(13) + 1.3 * parameter(13);
parameter(19) = parameter(19) + 1.8 * parameter(19);

parameter(1:13) = parameter(1:13) * b(34);
parameter(18:20) = parameter(18:20) * b(34);

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

PARAMETER = parameter;

NPOP = 24710;
CONS = 23.41;

rho = 0.2;

%tollini = [1;5];
tollini = 5.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A= [];
b = [];
Aeq = [];
beq = [];
lb = [0;2];
ub = [3;10];

options = optimset('TolFun', 1*10^-1, 'TolX', 1*10^-1, 'Display', 'iter', 'MaxFunEvals', 3000000);

tic
[tollexp, fval] = fminunc('optim1', tollini, options);
%toll = fmincon('optim2', tollini, A, b, Aeq, beq, lb, ub, 'nlcons', options);
%toll = fsolve(@optim2, tollini, options);
toc
stop
%% Calculate social welfare under equilibrium %%
tnew = [20.75;12.64];
tollnew = [1.98;8.35];

[res] = welfare_2step_new(NPOP, CONS, rho, parameter, tollnew, tnew, 8);

