function [tc] = timecost(tolle, tollf, te, tf, re, rf, scenario)

% scenario = 1: HOT
% scenario = 2: HOV
% scenario = 3: 2-route toll lanes

% We use the following labeling scheme
% 1: Solo_E
% 2: Solo_F

% 3: HOV2_E
% 4: HOV2_F

% 5: HOV3_E
% 6: HOV3_F

% Declare global variables
global CHOICE REPTITION NOBS NALT NFC NNC HIGHINC DISTANCE ...
       NERROR NFIXED NRANDOM RE RF KE KF SOLO HOV2 HOV3 SOLOE OCCUPANCE ...
       AFS PARAMETER NPOP CONS RHO; 

weights = NPOP / ( (NOBS/NALT) * REPTITION);
d2 = DISTANCE .* DISTANCE;
d3 = d2 .* DISTANCE;

%% Assigning time variables based on choice alternatives
timenew = zeros(NOBS,1) + te * RE + tf * RF;
relianew = zeros(NOBS,1) + re * RE + rf * RF;
if scenario == 1
	costnew = zeros(NOBS, 1) + tolle * SOLOE;
	masknew = [1;1;1;1;1;1];
elseif scenario == 2
	costnew = zeros(NOBS,1) + 10000 * SOLOE;
	masknew = [1;1;1;1;1;1];
else
	costnew = zeros(NOBS,1) + tolle * RE + tollf * RF;
	masknew = [1;1;1;1;1;1];
end;
costold = costold ./ OCCUPANCE;
costold_high = costold .* HIGHINC;
dtold = timeold .* DISTANCE;
dt2old = timeold .* d2;
dt3old = timeold .* d3;

costnew = costnew ./ OCCUPANCE;
costnew_high = costnew .* HIGHINC;
dtnew = timenew .* DISTANCE;
dt2new = timenew .* d2;
dt3new = timenew .* d3;

carpool = HOV2 + HOV3;
carpool_AFS = carpool .* AFS;

xfixedold = [costold,costold_high,dtold,dt2old,dt3old,reliaold,carpool,carpool_AFS,HOV3];
xrandomold = [timeold, HOV2, HOV3, reliaold];

xfixednew = [costnew,costnew_high,dtnew,dt2new,dt3new,relianew,carpool,carpool_AFS,HOV3];
xrandomnew = [timenew, HOV2, HOV3, relianew];

bfixed = PARAMETER(1:NFIXED);
brandom = PARAMETER( (NFIXED+1) : length(PARAMETER) );

fixold = xfixedold * bfixed;
fixnew = xfixednew * bfixed;

tao = zeros(NOBS/NALT,1);
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
    	indexold = fixold + brandom(1) * xold(:,1) + brandom(2) * (xold(:,2) + xold(:,3)) + ...
    	           bfixed(6) * brandom(3) * xold(:,4);
    	           
    	indexnew = fixnew + brandom(1) * xnew(:,1) + brandom(2) * (xnew(:,2) + xnew(:,3)) + ...
		   bfixed(6) * brandom(3) * xnew(:,4);
	
    	indexold = exp(indexold/RHO);    	
    	indexnew = exp(indexnew/RHO);
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tempold = zeros(NOBS/NALT, 1);
    	tempnew = zeros(NOBS/NALT, 1);
    	ini = 0;
    	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;
        	
        	if i == 1
        		tao(j,1) = PARAMETER(1)+PARAMETER(2)*HIGHINC(gstart);
               	end;

        	group = indexold(gstart:gend);
		group = group .* maskold;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(CONS + RHO*inclusive);
		tempold(j,1) = log(tscalar);
 		
        	group = indexnew(gstart:gend);
		group = group .* masknew;
	        denom = sum(group);
        	inclusive = log(denom);
        	tscalar = 1 + exp(CONS + RHO*inclusive);
        	tempnew(j,1) = log(tscalar);
        	
        	ini = ini + NALT;
        end;
	
	%tempev = tempnew - tempold;
	tempev = tempnew;
	tempev = (tempev) ./ (-1 * tao);
	if i == 1
		ev = tempev;
	else
		ev = [ev;tempev];
	end; 
end;

c75 = prctile(ev(:,1), 75);
c50 = prctile(ev(:,1), 50);
c25 = prctile(ev(:,1), 25);

w = sum(ev .* weights);


	
			
