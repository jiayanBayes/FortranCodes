function [f] = marginal(beta)
% This function calculates simulated likelihood function for mixed logit
% model

% Declare global variables
global REPTITION NUSERS NALT XFIXED XRANDOM YVEC NFC NNC ...
	CALPOLY SP WEIGHTS NERROR BROOK TIMES CHOICE RP ...
	SCALE_CAL SCALE_SP SCALE_BR; 

% Calculate the conditional probability for each individual and each
% time periods
bfixed = beta(1:NFC);
brandom = beta( (NFC+1) : (NFC+7) );
bscale = beta( (NFC+8) : length(beta));
fix = XFIXED * bfixed;

pn = zeros(NALT*NUSERS,1);
[r,c] = size(XRANDOM);

for i=1:REPTITION;
	x = zeros(r,c);
	count = 0;
	for j = 1:NNC;
    		g = count + i;
    		x(:,j) = XRANDOM(:,j) .* NERROR(:,g);
		count = count + REPTITION;
    	end;	 
   
    	index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) + brandom(3) * x(:,3) + brandom(4) * (x(:,4) + x(:,5))+ ...
    		brandom(5) * (RP .* x(:,6)) + brandom(6) * brandom(5) * (SP .* x(:,6)) + ...
    	  	bfixed(6) * brandom(7) * (RP .* x(:,7)) + ...
    		bfixed(22) * brandom(7) * (SP .* x(:,7));

    	%index = fix + brandom(1) * x(:,1) + brandom(2) * x(:,2) + brandom(3) * (x(:,3) + x(:,4))+ ...
    	%	brandom(4) * (RP .* x(:,5)) + brandom(5) * brandom(4) * (SP .* x(:,5));
    		
    	index = SP .* index + bscale(1) * (BROOK .* index) + ...
	    	bscale(2) * (CALPOLY .* index); 	    	
     
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% calculate conditional choice probality
    	ini = 0;
    	for j=1:(NUSERS);
        	gstart = ini + 1;
     		gend = ini + NALT;
     		group = index(gstart:gend);
     		group = exp(group);
     		denom = sum(group);
     		group = group ./ denom;
	
		if j == 1
			inter = group;
		else
			inter = [inter;group];
		end;

     	 	ini = ini + NALT;
   	 end;
	 clear index x;
   	 pn = pn + inter;
end;
f = pn / REPTITION;

	
			
