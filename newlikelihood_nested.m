function [f] = newlikelihood_nested(beta)
% This function calculates simulated likelihood function for mixed logit
% model

% Declare global variables
global REPTITION NUSERS NALT XFIXED XRANDOM YVEC NFC NNC ...
	CALPOLY SP WEIGHTS NERROR BROOK TIMES CHOICE RP ...
	SCALE_CAL SCALE_SP SCALE_BR; 

% Calculate the conditional probability for each individual and each
% time periods
bfixed = beta(1:NFC);
brandom = beta( (NFC+1) : (NFC+5) );
bscale = beta( (NFC+6) : (NFC+7) );
binc = beta( (NFC+8) : length(beta) );
fix = XFIXED * bfixed;

pn = zeros(NUSERS,1);
inter = zeros(NUSERS,1);
[r,c] = size(XRANDOM);

for i=1:REPTITION;
	x = zeros(r,c);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		x(:,j) = XRANDOM(:,j) .* NERROR(:,g);
		count = count + REPTITION;
    	end;	 
 
    	index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) +...
    		brandom(3) * (RP .* x(:,3)) + ...
    		brandom(4) * brandom(3) * (SP .* x(:,3)) + ...
    	  	bfixed(7)  * brandom(5) * (RP .* x(:,4)) + ...
    		bfixed(28) * brandom(5) * (SP .* x(:,4));
    		
    	index = SP .* index + bscale(1) * (BROOK .* index) + ...
	    	bscale(2) * (CALPOLY .* index); 	    	
     
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% calculate conditional choice probality
    	ini = 0;
    	for j=1:(NUSERS);
    		if TIMES(j) < NALT
   	 		% This case means that this person has only SP observations
    			inter(j) = 1;
    			for k = 1 : TIMES(j);
       				index(ini+k) = exp(-1*index(ini+k));
       				if YVEC(ini+k) == 1
       					inter(j) = inter(j)/(1+index(ini+k));
				else
					inter(j) = inter(j) * ...
					index(ini+k)/(1+index(ini+k));
       				end;
       			end;
   	 	elseif TIMES(j) == NALT
       			% This case means that this person has only RP observations
	      	  	gstart = ini + 1;
        			gend = ini + NALT;
        			group = index(gstart:gend);
				nest1 = group(1:3)/binc(1);
				nest2 = group(4:6)/binc(2);
				nest3 = group(7:9)/binc(3);
				snest1 = sum(nest1);
				snest2 = sum(nest2);
				snest3 = sum(nest3);
				logsum1 = log(snest1);
				logsum2 = log(snest2);
				logsum3 = log(snest3);
				denom = exp(binc(1) * logsum1) + exp(binc(2) * logsum2) + exp(binc(3) * logsum3);
				nest1 = nest1 ./ snest1;
				nest2 = nest2 ./ snest2;
				nest3 = nest3 ./ snest3;
				nest1 = nest1 * ( (exp(binc(1) * logsum1)) / denom );	
				nest2 = nest2 * ( (exp(binc(2) * logsum2)) / denom );	
				nest3 = nest3 * ( (exp(binc(3) * logsum3)) / denom );	

        			group = [nest1; nest2; nest3];
        		for k = gstart : gend;
		       		if YVEC(k) ==1 
            				inter(j) = group(k - ini);
         			 end;
      			end;
       		else
       			gstart = ini + 1;
       			gend = ini + NALT;
       			group = index(gstart:gend);
				nest1 = group(1:3)/binc(1);
				nest2 = group(4:6)/binc(2);
				nest3 = group(7:9)/binc(3);
				snest1 = sum(nest1);
				snest2 = sum(nest2);
				snest3 = sum(nest3);
				logsum1 = log(snest1);
				logsum2 = log(snest2);
				logsum3 = log(snest3);
				denom = exp(binc(1) * logsum1) + exp(binc(2) * logsum2) + exp(binc(3) * logsum3);
				nest1 = nest1 ./ snest1;
				nest2 = nest2 ./ snest2;
				nest3 = nest3 ./ snest3;
				nest1 = nest1 * ( (exp(binc(1) * logsum1)) / denom );	
				nest2 = nest2 * ( (exp(binc(2) * logsum2)) / denom );	
				nest3 = nest3 * ( (exp(binc(3) * logsum3)) / denom );	

     		   		group = [nest1; nest2; nest3];
       			for k = gstart : gend;
       				if YVEC(k) == 1
       					inter(j) = group(k - ini);
       				end;
       			end;
       			% For persons with both RP and SP observations, the first NALT
       			% observations are RP ones. Now we begin to calculate the
       			% conditional probability of SP observations from jth person 	
       			for k = 1 : (TIMES(j) - NALT);
       				index(gend+k) = exp(-1*index(gend+k));
       				if YVEC(gend+k) == 1
       					inter(j) = inter(j)/(1+index(gend+k));
				else
					inter(j) = inter(j) * ...
					index(gend+k)/(1+index(gend+k));
       				end;
       			end;
       		end;
      	 	ini = ini + TIMES(j);
   	 end;
	 clear index x;
   	 pn = pn + inter;
end;
pn = pn / REPTITION;
lnpn = log(pn);
lnpn = WEIGHTS .* lnpn; 
f = -1 * sum(lnpn);

	
			
