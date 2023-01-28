function [f] = func(beta)
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
bscale = beta( (NFC+6) : length(beta));
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
   
    	%index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) + brandom(3) * x(:,3) + brandom(4) * (x(:,4) + x(:,5))+ ...
    	%	brandom(5) * (RP .* x(:,6)) + brandom(6) * brandom(5) * (SP .* x(:,6)) + ...
    	%  	bfixed(6) * brandom(7) * (RP .* x(:,7)) + ...
    	%	bfixed(22) * brandom(7) * (SP .* x(:,7));
   
    	index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) + brandom(3) * (x(:,3) + x(:,4))+ ...
    		brandom(4) * (RP .* x(:,5)) + brandom(5) * brandom(4) * (SP .* x(:,5));
    		
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
        		group = exp(group);
        		denom = sum(group);
        		group = group ./ denom;
        		for k = gstart : gend;
		       		if YVEC(k) ==1 
            				inter(j) = group(k - ini);
         			 end;
      			end;
       		else
       			gstart = ini + 1;
       			gend = ini + NALT;
       			group = index(gstart:gend);
       			group = exp(group);
       			denom = sum(group);
       			group = group ./ denom;
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
f = WEIGHTS .* lnpn; 

	
			
