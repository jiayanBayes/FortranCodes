clear all;

% Define global variables
global REPTITION NOBS NFC NNC NALT NERROR NFIXED NRANDOM CHOICE ...
       HIGHINC DISTANCE NPOP CONS PARAMETER KE KF SOLO HOV2 HOV3 RE RF ...
       SOLOE OCCUPANCE AFS RHO FIXEDE FIXEDF STEPLENGTH DISCOUNT CONSTRUCTION ...
       WEIGHT BPOWER;  

%%%%%%%%%%% Set Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOMO = 0;

HALTON = 1;

REPTITION = 100;

DISCOUNT = 0.045;

CONSTRUCTION = 120000000; 

NOBS = 474;

NALT = 6;

NFC = 9;
NFIXED = 9;

NNC = 4;
NRANDOM = 4;

KE = 2000;
KF = 10000;

NPOP = 24300;
CONS = 23.41;
CONS = 42 - 4.6475 * (4-0.4);
RHO = 0.2;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load privateroad_data.txt;
data = privateroad_data;
clear privateroad_data;
CHOICE = data(:,2);
AFS = data(:,3);
DISTANCE = data(:,4);
HIGHINC = data(:,5);
MEDINC = data(:,6);
HIGHINC = HIGHINC + MEDINC;

SOLO    = zeros(NOBS,1);
HOV2    = zeros(NOBS,1);
HOV3    = zeros(NOBS,1);
RE 	= zeros(NOBS,1);
RF 	= zeros(NOBS,1);
OCCUPANCE = zeros(NOBS,1);
for i = 1:NOBS;
	if CHOICE(i) >= 1 & CHOICE(i) <= 2
		SOLO(i) = 1; 
		OCCUPANCE(i) = 1;
	end;
	if CHOICE(i) >= 3 & CHOICE(i) <=4
		HOV2(i) = 1;
		OCCUPANCE(i) = 2;
	end;
	if CHOICE(i) >= 5 & CHOICE(i) <= 6
		HOV3(i) = 1;
		OCCUPANCE(i) = 3; 
	end;
	if CHOICE(i) == 1 | CHOICE(i) == 3 | CHOICE(i) == 5
		RE(i) = 1; 
	end;
	if CHOICE(i) == 2 | CHOICE(i) == 4 | CHOICE(i) == 6
		RF(i) = 1; 
	end;
end;
SOLOE = RE .* SOLO;

load beta_revised_het.txt;
b = beta_revised_het;

beta = [b(1:6);b(13);b(15:16);b(28);b(30);b(33)];

% Rescale 
beta(1:9) = beta(1:9) * b(34);
beta(11)  = beta(11) * b(34);

%%%%%% Homogeneous Case %%%%%%%%%
if HOMO == 1

	NPOP = 27900;
	CONS = 11.82;
	CONS = 21.22 - 2.35* (4-0.4);
	tv1 = mean(DISTANCE);
	tv2 = mean(HIGHINC);
	tv3 = mean(AFS);
	
	for i = 1:length(DISTANCE);
		DISTANCE(i,1) = tv1;
		HIGHINC(i,1) = tv2;
		AFS(i,1) = tv3;
	end;
	beta(10) = 0;
	beta(11) = 0;
	beta(12) = 0;
end;
PARAMETER = beta;
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

%%%%%% No-Toll Equilibrium %%%%%%%%
STEPLENGTH = 0.03;
t1 = 0.00;
t2 = 9.84;
[junk, pre] = equilibrium(t1, t2, 3);

veupdate = sum(junk(:,1) + junk(:,3)/2 + junk(:,5)/3);
vfupdate = sum(junk(:,2) + junk(:,4)/2 + junk(:,6)/3);

neupdate = sum(junk(:,1) + junk(:,3)/1 + junk(:,5)/1);
nfupdate = sum(junk(:,2) + junk(:,4)/1 + junk(:,6)/1);

pfupdate = sum(junk(:,1) + junk(:,2));
solo = pfupdate / (sum(sum(junk)));
pfupdate = sum(junk(:,3) + junk(:,4));
shov2 = pfupdate / (sum(sum(junk)));
pfupdate = sum(junk(:,5)./1 + junk(:,6)./1);
shov3 = pfupdate / (sum(sum(junk))); 

tenotoll = supply(KE, veupdate)
tfnotoll = supply(KF, vfupdate)
renotoll = (tenotoll - 9.23) * 0.3785;
rfnotoll = (tfnotoll - 9.23) * 0.3785;

%[pa, pn1] = demand(0*1.1, 0*1.1, tenotoll*1.1, tfnotoll*1.1, renotoll*1.1, rfnotoll*1.1, 3);
%elas_fullcost = ( sum(sum(pn1)) - sum(sum(junk)) ) / sum(sum(junk))

%CONS = 21.22 - 2.35* 4 * 1.1;
%[pa, pn1] = demand(0, 0, tenotoll, tfnotoll, renotoll, rfnotoll, 3);
%elas_fuelcost = ( sum(sum(pn1)) - sum(sum(junk)) ) / sum(sum(junk))

tin = [t1;t2]
y = optim2(tin);
profit1 = (t1 * veupdate * 0.65 * 250 / DISCOUNT)
profit2 = (t2 * vfupdate * 0.65 * 250 / DISCOUNT) 

profit = profit1 + profit2  
stop
%%%Plot objective funcion %%%
%STEPLENGTH = 0.01
%x1 = linspace(0, 10, 5)'; 
%x2 = linspace(0, 10, 5)';
%counter = 0;
%dout = zeros(length(x1), length(x2));
%dppp = zeros(length(x1), length(x2));
%for i = 1:length(x1);
%	for j = 1:length(x2);
%		counter = counter + 1
%		tin = [x1(i);x2(j)];
%		y = optim2(tin);
%		dout(i,j) = y(1);
%		dppp(i,j) = y(2);
%	end;
%end;

%xi = linspace(0,10, 200);
%yi = linspace(0,10, 200)';
%w = interp2(x1, x2, dout, xi, yi);
%p = interp2(x1, x2, dppp, xi, yi);
%sortdata =  zeros(length(x1)*length(x2), 4);
%counter = 0;
%for i = 1:length(xi);
%	for j = 1:length(yi);
%		counter = counter + 1;
%		sortdata(counter,1) = xi(i);
%		sortdata(counter,2) = yi(j);
%		sortdata(counter,3) = w(i,j);
%		sortdata(counter,4) = p(i,j);
%	end;
%end;

%mesh(xi, yi, w); 
%hold on;
%mesh(xi, yi, p); 
%stop
%save surplus_heter.asc dout -ASCII;
%stop

%%% plot objective function in sequential game %%%
%STEPLENGTH = 0.03;
%points = linspace(1.47,12.92,100)';
%results = zeros(100,3);
%for i = 1:100;
%	i
%	tin = points(i,1);
%	results(i,1) = tin;
%	fval = optim1_sequential(tin);
%	results(i,2) = fval(1);
%	results(i,3) = fval(2);
%end;
%save welfare_leader_revision_privatebargaining.asc results -ASCII;
%stop

%%% OPTIMIZATION %%%
%% Solving response function of tolle with respect to tollf
%points = linspace(0, 10, 50)';
%results = zeros(50,2);
%for i = 1:50;
%	i
%	if i <= 10
%		tollini = 2;
%		STEPLENGTH = 0.015;
%	else
%		tollini = 2;
%		STEPLENGTH = 0.02;
%	end;
%	FIXEDF = points(i,1);
%	options = optimset('TolFun', 1*10^-3, 'TolX', 1*10^-3, 'Display', 'iter', 'MaxFunEvals', 3000000);
%	%[tolls, fval] = fminunc('optim1', tollini, options);
%	A= [];
%	b = [];
%	Aeq = [];
%	beq = [];
%	lb = [0];
%	ub = [20];
%	tolls = fmincon('optim1', tollini, A, b, Aeq, beq, lb, ub, 'nlcons', options);
%	results(i,1) = tolls;
%	results(i,2) = FIXEDF;
%	[tolls, FIXEDF]
%end;
%save response_operatorsolution_48.asc results -ASCII;
%stop

A= [];
b = [];
Aeq = [];
beq = [0];
%Aeq = [];
%beq = [];
lb = [0;0];
ub = [20;20];
STEPLENGTH = 0.03;
tollini = [1;10];
%tollini = 7;
options = optimset('TolFun', 1*10^-3, 'TolX', 1*10^-3, 'Display', 'iter', 'MaxFunEvals', 3000000);

points = linspace(1.50, 2.00, 10)';
results = zeros(10, 3);
%for i = 1:1;
%	BPOWER = points(i);
	%[tolls, fval] = fminunc('optim2', tollini, options);
	%[tolls, fval] = fmincon('optim2', tollini, A, b, Aeq, beq, lb, ub, 'nlcons', options);
	[tolls, fval] = fmincon('optim2', tollini, A, b, Aeq, beq, lb, ub);

%	results(i,1) = tolls(1);
%	results(i,2) = tolls(2);
%	results(i,3) = BPOWER;
%	results(i,:)
%	%toll = fsolve(@optim2, tollini, options);
%end;

