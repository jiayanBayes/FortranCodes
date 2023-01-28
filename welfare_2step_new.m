function [results] = welfare_2step_new(ntotal, cbar, rho, beta, tollnew, tnew, scenario)

% This function calculates the social welfare of different policy scenarios  
% The first three arguments: ntotal, cbar, and rho are three calibrated parameters.
% beta are is the vector of estimated parameters.
% tollnew and tnew are the vectors of time and toll of each alternative under proposed policy scenario.
% tollold and told are the vectors of time and toll of each alternative under base policy scenario.
% scenario is the proposed policy scenario.

% scenario = 0 means HOT1
% scenario = 1 means no-toll
% scenario = 2 means HOV
% scenario = 3 means HOT2 (only HOV3 are free)
% scenario = 4 means HOT3 (current scenario on SR91)
% scenario = 5 means general toll lanes
% scenario = 6 means two-route general differentiated toll
% scenario = 7 means two-route uniform
% scenario = 8 means two-route differentiated HOT
% scenario = 9 means two-route uniform HOT
 
%%%% Jia Yan 03/25/2005 %%%%

% Declare global variables
global CHOICE REPTITION NOBS NALT XMAT NFC NNC HIGHINC MEDINC DISTANCE...
       NERROR NFIXED NRANDOM IDIOSYNCRATIC OCCUPANCE EXPRESS SOLOEXP HOVEXP2 HOVEXP3...
       KE KF HOV2 HOV3 TRANS; 

weights = ntotal / ( (NOBS/NALT) * REPTITION);
bench = supply(KE, 0);
d2 = DISTANCE .* DISTANCE;
d3 = d2 .* DISTANCE;

% We use the following labeling scheme
% 1: Solo_notrans_91f
% 2: Solo_trans_91f
% 3: Solo_trans_91X

% 4: HOV2_notrans_91f
% 5: HOV2_trans_91f
% 6: HOV2_trans_91X

% 7: HOV3_notrans_91f
% 8: HOV3_trans_91f
% 9: HOV3_trans_91X

%% Assigning time variables based on choice alternatives
costold = zeros(NOBS,1);
timeold = zeros(NOBS,1) + 20.00 * EXPRESS + 20.00 * (ones(NOBS,1) - EXPRESS);

maskold = [1;0;0;1;0;0;1;0;0];
reliaold = 0.3785 * (timeold - bench * ones(NOBS,1));

timenew = zeros(NOBS,1) + tnew(2) * EXPRESS + tnew(1) * (ones(NOBS,1) - EXPRESS);
rnew = 0.3785 * (tnew - bench * ones(2,1));
relianew = zeros(NOBS,1) + rnew(2) * EXPRESS + rnew(1) * (ones(NOBS,1) - EXPRESS);


if scenario == 0
	costnew = zeros(NOBS, 1) + tollnew(2) * SOLOEXP;
	masknew = [1;1;1;1;1;1;1;1;1];
elseif scenario == 1
	costnew = zeros(NOBS,1);
	masknew = [1;0;0;1;0;0;1;0;0];
elseif scenario == 2
	costnew = zeros(NOBS,1) + 10000 * SOLOEXP;;
	masknew = [1;0;0;1;0;1;1;0;1];
elseif scenario == 3
	costnew = zeros(NOBS,1) + tollnew(2) * (SOLOEXP + HOVEXP2) + 0 * tollnew(2) * HOVEXP3;
	masknew = [1;1;1;1;1;1;1;1;1];
elseif scenario == 4
	costnew = zeros(NOBS,1) + tollnew(2) * (SOLOEXP + HOVEXP2) + 0.5 * tollnew(2) * HOVEXP3;
	masknew = [1;1;1;1;1;1;1;1;1];
elseif scenario == 5
	costnew = zeros(NOBS,1) + tollnew(2) * (SOLOEXP + HOVEXP2 + HOVEXP3);
	masknew = [1;1;1;1;1;1;1;1;1];
elseif scenario == 6
	costnew = zeros(NOBS,1) + tollnew(2) * (SOLOEXP + HOVEXP2 + HOVEXP3) + ...
		  tollnew(1) * (ones(NOBS,1) - EXPRESS);
	masknew = [1;0;1;1;0;1;1;0;1];
elseif scenario == 7
	costnew = zeros(NOBS,1) + tollnew(2) * (SOLOEXP + HOVEXP2 + HOVEXP3) + ...
		  tollnew(1) * (ones(NOBS,1) - EXPRESS);
	masknew = [1;0;0;1;0;0;1;0;0];
elseif scenario == 8
	costnew = zeros(NOBS,1) + tollnew(2) * SOLOEXP + ...
		  tollnew(1) * (ones(NOBS,1) - SOLOEXP - HOV2 - HOV3);
	masknew = [1;0;1;1;0;1;1;0;1];
else
	costnew = zeros(NOBS,1) + tollnew(2) * SOLOEXP + ...
		  tollnew(1) * (ones(NOBS,1) - SOLOEXP - HOV2 - HOV3);
	masknew = [1;0;0;1;0;0;1;0;0];
end;
costold = costold ./ OCCUPANCE;
costold_high = costold .* HIGHINC;
costold_med = costold .* MEDINC;

dtold = timeold .* DISTANCE;
dt2old = timeold .* d2;
dt3old = timeold .* d3;

costnew = costnew ./ OCCUPANCE;
costnew_high = costnew .* HIGHINC;
costnew_med = costnew .* MEDINC;

dtnew = timenew .* DISTANCE;
dt2new = timenew .* d2;
dt3new = timenew .* d3;

xfixedold = [costold,costold_high,dtold,dt2old,dt3old,reliaold,XMAT];
xrandomold = [timeold, TRANS, HOV2, HOV3, EXPRESS, reliaold];

xfixednew = [costnew,costnew_high,dtnew,dt2new,dt3new,relianew,XMAT];
xrandomnew = [timenew, TRANS, HOV2, HOV3, EXPRESS, relianew];

bfixed = beta(1:NFIXED);
bfixed_ntrans = bfixed;
brandom = beta( (NFIXED+1) : length(beta) );

%% to get the base case %%
bfixed_ntrans(7:9) = 0;
bfixed_ntrans(10) = bfixed_ntrans(10) - 0.1552; % Base case
bfixed_ntrans(14:16) = 0;

fixold_ntrans = xfixedold * bfixed_ntrans;

%% Do the adjustments based on Ken's suggestions %%
if scenario == 6 | scenario == 8
	xfixedold(:,7) = 1;
	xfixednew(:,7) = 1;

	ini = 0;
	%% reset the interactions with transponder and express lane dummies %%
	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;

        	group = xfixedold(gstart:gend, 8:9);
       		xfixedold(gstart:gend,8) = group(3,1);
       		xfixedold(gstart:gend,9) = group(3,2);
        	
        	group = xfixednew(gstart:gend, 8:9);
       		xfixednew(gstart:gend,8) = group(3,1);
       		xfixednew(gstart:gend,9) = group(3,2);
        	
		ini = ini + NALT;
	end;
	xrandomold(:,2) = 1;
	xrandomnew(:,2) = 1;
end;

if scenario == 2 
	bfixed(7:9) = 0;
	brandom(2) = 0;
end;

fixold = xfixedold * bfixed;
fixnew = xfixednew * bfixed;

tao = zeros(NOBS/NALT,1);
hincome = zeros(NOBS/NALT, 1);
for i=1:REPTITION;
	xold = zeros(NOBS, NNC);
	xnew = zeros(NOBS,NNC);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		xold(:,j) = xrandomold(:,j) .* NERROR(:,g);    		
    		xnew(:,j) = xrandomnew(:,j) .* NERROR(:,g);    		
		count = count + REPTITION;
    	end;	 
    	indexold_ntrans = fixold_ntrans + brandom(1) * xold(:,1) + 0 * brandom(2) * xold(:,2) + ...
    	           brandom(3) * (xold(:,3) + xold(:,4)) + ...
    	           0 * brandom(4) * xold(:,5) + bfixed(6) * brandom(5) * xold(:,6);

    	indexold = fixold + brandom(1) * xold(:,1) + brandom(2) * xold(:,2) + ...
    	           brandom(3) * (xold(:,3) + xold(:,4)) + ...
    	           brandom(4) * xold(:,5) + bfixed(6) * brandom(5) * xold(:,6);

    	indexnew = fixnew + brandom(1) * xnew(:,1) + brandom(2) * xnew(:,2) + ...
    	           brandom(3) * (xnew(:,3) + xnew(:,4)) + ...
    	           brandom(4) * xnew(:,5) + bfixed(6) * brandom(5) * xnew(:,6);
	
    	indexold_ntrans = exp(indexold_ntrans/rho);    	
    	indexold = exp(indexold/rho);    	
    	indexnew = exp(indexnew/rho);

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tempold = zeros(NOBS/NALT, 1);
    	tempnew1 = zeros(NOBS/NALT, 1);
    	tempnew2 = zeros(NOBS/NALT, 1);
    	tempnew3 = zeros(NOBS/NALT, 1);

    	ini = 0;
    	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;
        	
        	if i == 1
        		tao(j,1) = beta(1)+beta(2)*HIGHINC(gstart);
        		hincome(j,1) = HIGHINC(gstart);
               	end;

		% the base
        	group = indexold_ntrans(gstart:gend);
		group = group .* maskold;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(cbar + rho*inclusive);
        	
        	%% Elastic case %%
        	tempold(j,1) = log(tscalar);
		
		%% Inelastic Case %%
		%tempold(j,1) = inclusive;
		
		% first step: adding new alternatives
        	group = indexold_ntrans(gstart:gend);
		group = group .* masknew;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(cbar + rho*inclusive);

		% Elastic Case %
        	tempnew1(j,1) = log(tscalar);
		
		%% Inelastic Case %%
		%tempnew1(j,1) = inclusive;

		% second step: adding preference for transponder and lane
        	group = indexold(gstart:gend);
		group = group .* masknew;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(cbar + rho*inclusive);
        	tempnew2(j,1) = log(tscalar);
        	%tempnew2(j,1) = inclusive;
        	
        	% third step: time and cost change
        	group = indexnew(gstart:gend);
		group = group .* masknew;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(cbar + rho*inclusive);
        	tempnew3(j,1) = log(tscalar);
        	%tempnew3(j,1) = inclusive;
        	
        	ini = ini + NALT;
        end;
	
	tempev1 = tempnew1 - tempold;
	tempev2 = tempnew2 - tempnew1;
	tempev3 = tempnew3 - tempnew2;
	tempev4 = tempnew3 - tempold;
	
	tempev1 = (tempev1) ./ (-1 * tao);
	tempev2 = (tempev2) ./ (-1 * tao);
	tempev3 = (tempev3) ./ (-1 * tao);
	tempev4 = (tempev4) ./ (-1 * tao);
	

	if i == 1
		ev1 = tempev1;
		ev2 = tempev2;
		ev3 = tempev3;
		ev4 = tempev4;
		inc = hincome;
	else
		ev1 = [ev1;tempev1];
		ev2 = [ev2;tempev2];
		ev3 = [ev3;tempev3];
		ev4 = [ev4;tempev4];
		inc = [inc; hincome];
	end;
end;

c75 = prctile(ev4(:,1), 75);
c50 = prctile(ev4(:,1), 50);
c25 = prctile(ev4(:,1), 25);

wstep1 = sum(ev1 .* weights);
wstep2 = sum(ev2 .* weights); 
wstep3 = sum(ev3 .* weights); 
wstep4 = sum(ev4 .* weights);

% distribution by income group
dim1 = floor(sum(inc));
dim2 = floor(sum(ones(length(inc),1) - inc));
ev5 = zeros(dim1, 1);
ev6 = zeros(dim2, 1);

num1 = 0;
num2 = 0;
for i = 1:length(ev4);
	if inc(i,1) == 1
		num1 = num1 + 1;
		ev5(num1,1) = ev4(i,1);
	else
		num2 = num2 + 1;
		ev6(num2,1) = ev4(i,1);
	end;
end;

hc75 = prctile(ev5, 75);
hc50 = prctile(ev5, 50);
hc25 = prctile(ev5, 25);

lc75 = prctile(ev6, 75);
lc50 = prctile(ev6, 50);
lc25 = prctile(ev6, 25);

results = [wstep1; wstep2; wstep3; wstep4; c75; c50; c25; hc75; hc50; hc25; lc75; lc50; lc25];




	
			
