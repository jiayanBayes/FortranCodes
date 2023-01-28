function [f] = mixed_logit_rp(beta)
% This function calculates simulated likelihood function for mixed logit
% model
% Declare global variables
global REPTITION NOBS NALT XMAT YVEC NFC IDFC NNC IDNC CHOICE...
	CALPOLY BROOK WEIGHTS HRANDN HRANDN_SCALE NERROR HALTON; 

% Calculate the conditional probability for each individual and each
% time periods
pn = zeros(NOBS/NALT, 1);
inter = zeros(NOBS/NALT, 1);
for i=1:REPTITION;
    x=zeros(NOBS,(NFC+2*NNC));
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
      %x1 = x(:,1:16);
      %x2 = x(:,18);
      %nx = [x1,x2];
     index = x * beta(1:16);
     
     for j = 1:NOBS;
     	if CALPOLY(j) == 1
     		index(j) = beta(17) * index(j);
     	end;			
     end;
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate conditional choice probality
    index=exp(index);
    ini = 0;
    for j=1:(NOBS/NALT);
        gstart = ini + 1;
        gend = ini + NALT;
        group = index(gstart:gend);
        denom = sum(group);
        group = group ./ denom;
        for k = gstart : gend;
          if YVEC(k) ==1 
            inter(j) = group(k - ini);
          end;
       end;
       ini = ini + NALT;
     end;
    clear index x;
    pn = pn + inter;
end;

pn = pn / REPTITION;
lnpn = log(pn);
lnpn = WEIGHTS .* lnpn; 
f = -1 * sum(lnpn);

	
			
