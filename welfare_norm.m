function [wel1, wel2, ev4] = welfare_norm(ntotal, cbar, rho, beta, tollnew, tollold, tnew, told, scenario)

% This function calculates the social welfare of different policy scenarios  
% The first three arguments: ntotal, cbar, and rho are three calibrated parameters.
% beta are is the vector of estimated parameters.
% tollnew and tnew are the vectors of time and toll of each alternative under proposed policy scenario.
% tollold and told are the vectors of time and toll of each alternative under base policy scenario.
% scenario is the proposed policy scenario.

%% We take HOT1 (all HOVs are free) as the base case %%
% scenario = 0: HOT1 and By-pass lanes
% scenario = 1: HOT1 and No-toll
% scenario = 2: HOT1 and HOV
% scenario = 3: HOT1 and HOT2 (only HOV3 are free)
% scenario = 4: HOT1 and HOT3 (HOV3 pay half)
% scenario = 5: HOT1 and genreal toll lanes
% scenario = 6: HOT1 and two-route general differentiated toll
% scenario = 7: HOT1 and two-route uniform toll or the first-best toll
% scenario = 8: HOT1 and two-route differentiated HOT
% scenario = 9: HOT1 and two-route uniform HOT

%%%% Jia Yan 03/25/2005 %%%%

% Declare global variables
global CHOICE REPTITION NOBS NALT XMAT NFC NNC HIGHINC MEDINC DISTANCE...
       NERROR NFIXED NRANDOM IDIOSYNCRATIC OCCUPANCE EXPRESS SOLOEXP HOVEXP2 HOVEXP3...
       KE KF HOV2 HOV3 TRANS; 

weights = ntotal / ( (NOBS/NALT) * REPTITION);
bench = supply(KE, 0);
d2 = DISTANCE .* DISTANCE;
d3 = d2 .* DISTANCE;

% We use the following labeling scheme for choice alternatives
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
timenew = zeros(NOBS,1) + tnew(2) * EXPRESS + tnew(1) * (ones(NOBS,1) - EXPRESS);
timeold = zeros(NOBS,1) + told(2) * EXPRESS + told(1) * (ones(NOBS,1) - EXPRESS);

rnew = 0.3785 * (tnew - bench * ones(2,1));
rold = 0.3785 * (told - bench * ones(2,1)); 
relianew = zeros(NOBS,1) + rnew(2) * EXPRESS + rnew(1) * (ones(NOBS,1) - EXPRESS);
reliaold = zeros(NOBS,1) + rold(2) * EXPRESS + rold(1) * (ones(NOBS,1) - EXPRESS);

%% Asisgning tolls and others based on alternatives and policy scenarios %%

%% The base case is HOT1 %% 
costold = zeros(NOBS,1) + tollold(2) * SOLOEXP + 0 * tollold(2) * (HOVEXP2 + HOVEXP3); 
maskold = [1;0;1;1;0;1;1;0;1];
if scenario == 0
	costnew = zeros(NOBS,1);
	masknew = maskold
elseif scenario == 1
	costnew = zeros(NOBS,1);
	masknew = [1;0;0;1;0;0;1;0;0];
elseif scenario == 2
	costnew = zeros(NOBS,1);
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

%% we assume that the carpoolers share their money costs %% 
costnew = costnew ./ OCCUPANCE;
costold = costold ./ OCCUPANCE;

%% Interact cost and time with users' profiles %%
costold_high = costold .* HIGHINC;
costnew_high = costnew .* HIGHINC;
costold_med =  costold .* MEDINC;
costnew_med =  costnew .* MEDINC;

dtold =  timeold .* DISTANCE;
dt2old = timeold .* d2;
dt3old = timeold .* d3;
dtnew =  timenew .* DISTANCE;
dt2new = timenew .* d2;
dt3new = timenew .* d3;


%% To obtain the X matrix %%
xfixedold = [costold,costold_high,costold_med,dtold,dt2old,dt3old,reliaold,XMAT];
xfixednew = [costnew,costnew_high,costnew_med,dtnew,dt2new,dt3new,relianew,XMAT];

xrandomold = [timeold,TRANS,HOV2,HOV3,EXPRESS,reliaold];
xrandomnew = [timenew,TRANS,HOV2,HOV3,EXPRESS,relianew];

%% To get the non-random part of XB %%
bfixed =  beta(1:NFIXED);
brandom = beta( (NFIXED+1) : length(beta) );
bfixed_ntrans = bfixed;
bfixed_ntrans(8:13) = 0; % reset the transponder dummy as zero

fixold =  xfixedold * bfixed;
fixold_ntrans = xfixedold * bfixed_ntrans;
if scenario <= 2 | scenario >= 6
	bfixed_new = bfixed_ntrans;
	bfixed_new(14) = bfixed_new(14) - 0.795;
else
	bfixed_new = bfixed;
end;

fixold =  xfixedold * bfixed;
fixold_ntrans = xfixedold * bfixed_ntrans;
fixnew =  xfixednew * bfixed_new;

tao = zeros(NOBS/NALT,1);
for i=1:REPTITION;
	xold = zeros(NOBS,NNC);
	xnew = zeros(NOBS,NNC);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		xold(:,j) = xrandomold(:,j) .* NERROR(:,g);
    		xnew(:,j) = xrandomnew(:,j) .* NERROR(:,g);    		
		count = count + REPTITION;
    	end;	 
   
    	indexold = fixold + brandom(1) * xold(:,1) + brandom(2) * xold(:,2) + ...
    	           brandom(3) * (xold(:,3) + xold(:,4)) + ...
    	           brandom(4) * xold(:,5) + bfixed(7) * brandom(5) * xold(:,6);
    	
    	indexold_ntrans = fixold_ntrans + brandom(1) * xold(:,1) + brandom(2) * xold(:,2) + ...
    			  brandom(3) * (xold(:,3) + xold(:,4)) + ...
    			  brandom(4) * xold(:,5) + bfixed(7) * brandom(5) * xold(:,6);
    	
    	indexnew = fixnew + brandom(1) * xnew(:,1) + brandom(2) * xnew(:,2) + ...
    	           brandom(3) * (xnew(:,3) + xnew(:,4)) + ...
    	           brandom(4) * xnew(:,5) + bfixed(7) * brandom(5) * xnew(:,6);
    	    	
     	indexold = exp(indexold);
     	indexold_ntrans = exp(indexold_ntrans);
    	indexnew = exp(indexnew);
    	
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	tempold1 = zeros(NOBS/NALT, 1);
    	tempnew1 = zeros(NOBS/NALT, 1);
    	tempnew2 = zeros(NOBS/NALT, 1);
    	tempnew3 = zeros(NOBS/NALT, 1);
    	tempnew4 = zeros(NOBS/NALT, 1);
    	ini = 0;
    	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;
        	
        	%% To obtain the marginal utility of income for each person %%
        	if i == 1
        		tao(j,1) = beta(1)+beta(2)*HIGHINC(gstart)+beta(3)*MEDINC(gstart);
               	end;
               	
		%%%%% The 1st Step of welfare calculation: eliminating transponder choices %%%% 
		group =  indexold(gstart:gend);
               	denom = sum(group);
               	inclusive = log(denom);
        	tscalar = 1 + exp(cbar+rho*inclusive);
        	%tempold1(j,1) = log(tscalar);
        	tempold1(j,1) = inclusive;

		group =  indexold(gstart:gend);
		group = group .* maskold;
               	denom = sum(group);
               	inclusive = log(denom);
        	tscalar = 1 + exp(cbar+rho*inclusive);
        	%tempnew1(j,1) = log(tscalar);
        	tempnew1(j,1) = inclusive;
 		
		%%%%% The 2nd Step of welfare calculation: further eliminating transponder preferences %%%% 
 		group =  indexold_ntrans(gstart:gend);
 		group = group .* maskold;
               	denom = sum(group);
               	inclusive = log(denom);
        	tscalar = 1 + exp(cbar+rho*inclusive);
        	%tempnew2(j,1) = log(tscalar);
        	tempnew2(j,1) = inclusive;
        	
		%% The 3rd step of welfare calculation: To obtain the "log-sum" of the proposed case %%
        	group = indexnew(gstart:gend);
        	group = group .* masknew;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(cbar + rho*inclusive);
        	%tempnew3(j,1) = log(tscalar);
        	tempnew3(j,1) = inclusive;

        	ini = ini + NALT;
        end;
	
	tempev1 = (tempnew1 - tempold1);
	tempev1 = tempev1 ./ (-1*tao);

	tempev2 = (tempnew2 - tempnew1);
	tempev2 = tempev2 ./ (-1*tao);

	tempev3 = (tempnew3 - tempnew2);
	tempev3 = tempev3 ./ (-1*tao);
	
	tempev4 = tempev3 + tempev2 + tempev1;
	
	if i == 1
		ev1 = tempev1;
		ev2 = tempev2;
		ev3 = tempev3;
		ev4 = tempev4;
	else
		ev1 = [ev1;tempev1];
		ev2 = [ev2;tempev2];
		ev3 = [ev3;tempev3];
		ev4 = [ev4;tempev4];
	end; 
end;
wel1 = sum(ev3 .* weights);
wel2 = sum(ev4 .* weights);
