function [pa, pn] = demand(tolle, tollf, te, tf, re, rf, scenario)
% This function calculates the demand for different routes, and modes
%% We consider different scenarios %%
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
       NERROR NFIXED NRANDOM RE RF KE KF SOLO HOV2 HOV3 SOLOE OCCUPANCE...
       AFS NPOP CONS RHO PARAMETER; 

weights = NPOP / ( (NOBS/NALT) * REPTITION);

d2 = DISTANCE .* DISTANCE;
d3 = d2 .* DISTANCE;

%% Assigning time variables based on choice alternatives
time  = zeros(NOBS,1) + te * RE + tf * RF;
relia = zeros(NOBS,1) + re * RE + rf * RF;

%% Setting tolls and others based on different policy scenarios%%
if scenario == 1
	cost = zeros(NOBS,1) + tolle * SOLOE;
elseif scenario == 2
	cost = zeros(NOBS,1) + 1000000 * SOLOE;
else
	cost = zeros(NOBS,1) + tolle * RE + tollf * RF;
end;
cost = cost ./ OCCUPANCE;

%% Construct cost and time variables in the model
cost_high = cost .* HIGHINC;
dt = time .* DISTANCE;
dt2 = time .* d2;
dt3 = time .* d3;
carpool = HOV2 + HOV3;
carpool_AFS = carpool .* AFS;

xfixed = [cost, cost_high, dt, dt2, dt3, relia, carpool, carpool_AFS, HOV3];
xrandom = [time, HOV2, HOV3, relia];

bfixed = PARAMETER(1:NFIXED);
brandom = PARAMETER( (NFIXED+1) : length(PARAMETER) );
fix = xfixed * bfixed;

temp1 = zeros(NOBS/NALT,1);
temp2 = zeros(NOBS/NALT,NALT);
junk = zeros(NOBS/NALT,2);
for i=1:REPTITION;
	x = zeros(NOBS,NNC);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		x(:,j) = xrandom(:,j) .* NERROR(:,g);
		count = count + REPTITION;
    	end;	 
   
    	index = fix + brandom(1) * x(:,1)+ brandom(2) * (x(:,2) + x(:,3)) + ...
   	        bfixed(6) * brandom(3) * x(:,4);
   	        
    	index = exp(index/RHO);

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	ini = 0;
    	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;

       		group = index(gstart:gend);
	       	denom = sum(group);
       		inclusive = log(denom);
       		group = group ./ denom;
        		
       		temp1(j,1) = exp(CONS+RHO*inclusive)/(1+exp(CONS+RHO*inclusive));
       	       	       	
       		for k = 1 : NALT;
			temp2(j,k) = group(k) * temp1(j,1);
    		end;
      		
               	ini = ini + NALT;
	end;
	if i == 1
		pa = temp1;
		pn = temp2;
	else
		pa = [pa;temp1];
		pn = [pn;temp2];
	end;
end;
% pa is the demand for highway travel; pn is the demand for each alternative. 
pa = weights * pa;
pn = weights * pn;




	
			
