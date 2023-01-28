function [pa,pn] = welfare_1step(ntotal, cbar, rho, beta, tolle, tollf, te, tf, re, rf, tup, scenario)
% This function calculates the demand for different routes, and modes
%% We consider different scenarios %%
% scenario = 1 means no-toll case;
% scenario = 2 means HOV case;
% scenario = 3 means HOT case and HOVs are free;
% scenario = 4 means HOT case and HOV3 pay only half;
% scenario = 5 means general Toll lanes; 

% Declare global variables
global CHOICE REPTITION NOBS NALT XMAT NFC NNC HIGHINC MEDINC DISTANCE ...
       NERROR NFIXED NRANDOM IDIOSYNCRATIC EXPRESS OCCUPANCE SOLOEXP HOVEXP2 HOVEXP3; 

weights = ntotal / ( (NOBS/NALT) * REPTITION);

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
time  = zeros(NOBS,1) + te * EXPRESS + tf * (ones(NOBS,1) - EXPRESS);
relia = zeros(NOBS,1) + re * EXPRESS + rf * (ones(NOBS,1) - EXPRESS);

%% Setting tolls and others based on different policy scenarios%%
if scenario < 3
	cost = zeros(NOBS,1);
elseif scenario == 3
	cost = zeros(NOBS,1) + tolle * SOLOEXP;
elseif scenario == 4
	cost = zeros(NOBS,1) + tolle * (SOLOEXP + HOVEXP2) + 0 * tolle * HOVEXP3;
else
	cost = zeros(NOBS,1) + tolle * (SOLOEXP + HOVEXP2 + HOVEXP3);
end;
cost = cost ./ OCCUPANCE;

%% Construct cost and time variables in the model
cost_high = cost .* HIGHINC;
cost_med = cost .* MEDINC;
dt = time .* DISTANCE;
dt2 = time .* d2;
dt3 = time .* d3;


xfixed = [cost, cost_high, cost_med, dt, dt2, dt3, relia, XMAT];
xrandom = [time,EXPRESS,relia];

bfixed = beta(1:NFIXED);
brandom = beta( (NFIXED+1) : (NFIXED+NRANDOM) );

fix = xfixed * bfixed;

%% for scenario 2 (HOV case), we force that its is impossible for solo drivers
%% to use express lanes 
if scenario == 2 
	fix = fix - 10000 * SOLOEXP;
end;

temp1 = zeros(NOBS/NALT,1);
temp2 = zeros(NOBS/NALT,2);

bmark = supply(4000, 0);
rup = 0.3785 * (tup - bmark);
for i=1:REPTITION;
	x = zeros(NOBS,NNC);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		x(:,j) = xrandom(:,j) .* NERROR(:,g);
		count = count + REPTITION;
    	end;	 
   
    	index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) + bfixed(7) * brandom(3) * x(:,3);
    	
    	vt1 = bfixed(1) * ones(length(HIGHINC),1) + bfixed(2) * HIGHINC + bfixed(3) * MEDINC;
    	vt2 = bfixed(4) * DISTANCE + bfixed(5) * d2 + bfixed(6) * d3 + brandom(1) * NERROR(:,i);
    	vt3 = bfixed(7) * ones(length(HIGHINC),1) + bfixed(7) * brandom(3) * NERROR(:, (REPTITION + i));
    	vtime = (vt2 * 60) ./ vt1;
    	vr = (vt3 * 60) ./ vt1;
       	%index = index/rho;
    	  	
    	index = exp(index);
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	ini = 0;
    	for j = 1:(NOBS/NALT); 
    		gstart = ini + 1;
        	gend = ini + NALT;

       		group = index(gstart:gend);
	       	denom = sum(group);
       		inclusive = log(denom);
       		group = group ./ denom;
        	
       		temp1(j,1) = exp(cbar+rho*inclusive)/(1+exp(cbar+rho*inclusive));
       	       	       	
       		for k = 1 : NALT;
			group(k) = group(k)* temp1(j,1);
       		end;
       		
       		pe = group(3) + group(6) + group(9);
       		pf = group(1) + group(2) + group(4) + group(5) + group(7) + group(8);
       		
       		temp2(j, 1) = pe * vtime(gstart,1) * (tup - te) + pf * vtime(gstart,1) * (tup - tf);
       		temp2(j, 2) = pe * vr(gstart,1) * (rup - re) + pf * vr(gstart,1) * (rup - rf);
       		
               	ini = ini + NALT;
	end;
	
	if i == 1
		wt = temp2;
	else
		wt = [wt;temp2];
	end;
end;
pa = weights * pa;
pn = weights * pn;




	
			
