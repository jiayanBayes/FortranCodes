function [f] = likelihood_func(beta)
% This function calculates simulated likelihood function for mixed logit
% model
% Declare global variables
global REPTITION NUSERS NALT XMAT YVEC NFC IDFC NNC IDNC...
	CALPOLY SP WEIGHTS HALTON HRANDN HRANDN_SCALE_CAL...
	HRANDN_SCALE_SP HRANDN_SCALE_BR NERROR BROOK TIMES CHOICE RP; 

% Calculate the conditional probability for each individual and each
% time periods
pn = zeros(NUSERS, 1);
inter = zeros(NUSERS, 1);
for i=1:REPTITION;
    x=zeros(sum(TIMES),(NFC+2*NNC));
	% Add fixed variables	
	if NFC>0
		for j=1:NFC;
			x(:,j)=XMAT(:,IDFC(j));
		end;
	else
		x=x;
	end;
	% Add normal distributed variables
	if NNC>0
       		if HALTON == 1
        	    	iter1=1;
	    		for j=1:NNC;
				nstart = (j-1)*REPTITION + i;
				x(:,NFC+iter1) = XMAT(:,IDNC(j));
				x(:,NFC+iter1+1) = HRANDN(:,nstart).*XMAT(:,IDNC(j));
				iter1=iter1+2;
	    		end;	
        	else
            		iter = 1;
            		for j=1:NNC;
            			nstart = (j-1)*REPTITION + i;
				x(:,NFC+iter) = XMAT(:,IDNC(j));
				x(:,NFC+iter+1) = NERROR(:,nstart) .* XMAT(:,IDNC(j));
				iter = iter + 2;
	    		end;
        	end;
      else
		x=x;
      end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This part will reorganize X matrix based on specific model
    x1 = x(:,1:36);
    x2 = x(:,38);
    x3 = x(:,40);
    x4 = x(:,42);
    x5 = x(:,44);
    x6 = x(:,46);
    x7 = x(:,48);
    x8 = x(:,50);
    x9 = x(:,52);
    x10 = x(:,54); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = [x1, x2, x3, x4, x5, x6];
    index = nx * beta(1:41) + ...
    		BROOK .* (x7 + x8 + x9) + beta(42) * (SP.*(x7 + x8 + x9)) + ...
    		beta(5) * beta(43)* (RP .* x10) + ...
    		beta(30) * beta(43) * (SP .* x10) +...
    		BROOK .* HRANDN_SCALE_BR(:,i) + ...
    		beta(44) * (SP .* HRANDN_SCALE_SP(:,i)) + ...
    		beta(45) * (CALPOLY .* HRANDN_SCALE_CAL(:,i));
  
    %index = x * beta(1:34);
    
    %index = beta(35) * (SP .* index) + beta(36) * (CALPOLY .* index) + ...
    		BROOK .* index;
  
    %index = x * beta(1:34) + beta(35) * (SP .* HRANDN_SCALE_SP(:,i)) + ...
    		beta(36) * (CALPOLY .* HRANDN_SCALE_CAL(:,i));
    
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
%f = -1 * sum(lnpn);

	
			
