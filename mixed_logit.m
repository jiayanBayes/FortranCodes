function [f] = mixed_logit(beta)
% This function calculates simulated likelihood function for mixed logit
% model
% Declare global variables
global REPTITION NOBS NALT NG XMAT YVEC TIMES NFC IDFC NNC IDNC...
	REVEAL STATED CALPOLY BROOK WEIGHTS; 

% Calculate the conditional probability for each individual and each
% time periods
pnt = zeros(length(NALT), 1);
pn = zeros(NG, 1);
inter = zeros(NG, 1);
seed = 13457;
for i=1:REPTITION;
    seed = seed + i;
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
        randn('state', seed);
        draw = randn(NG, NNC);
        ini=0;
        NERROR = zeros(NOBS, NNC);
        for j=1:NG;
        	gstart=ini+1;
          	gend=ini+TIMES(j);
            temp0=ones(TIMES(j),1);
            temp1=kron(draw(j,:),temp0);
            NERROR(gstart:gend,:)=temp1;
            ini=ini+TIMES(j);
        end;
       
        iter = 1;
        for j=1:NNC;
			x(:,NFC+iter) = XMAT(:,IDNC(j));
			x(:,NFC+iter+1) = NERROR(:,j) .* XMAT(:,IDNC(j));
			iter = iter + 2;
		end;	
    else
		x=x;
	end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This part will reorganize X matrix based on specific model
    
 	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	index = nx * beta;
	for j=1:NOBS
		if REVEAL(j)==1
			index(j) = index(j) + nx(j,22);
		else
			index(j) = index(j) + beta(22) * nx(j,22);
            index(j) = scale_sp * index(j);
		end;
            
        if CALPOLY(j) == 1
            index(j) = scale_cal * index(j);
        end;
        
	end;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	% calculate conditional choice probality
    index=exp(index);
    ini = 0;
    for j=1:length(NALT);
        gstart = ini + 1;
        gend = ini + NALT(j);
        group = index(gstart:gend);
        denom = sum(group);
        group = group ./ denom;
        for k = gstart : gend;
          if YVEC(k) ==1 
            pnt(j) = group(k - ini);
          end;
       end;
       ini = ini + NALT(j);
     end;
        
	clear index x;
    
    ini = 0;
    for j = 1: NG;
        gstart = ini + 1;
        gend = ini + TIMES(j);
        group = pnt(gstart:gend);
        inter(j) = prod(group);
        ini = ini + TIMES(j)
    end;
    pn = pn + inter;
end;
pn = pn ./ REPTITION;

lnpn = log(pn);
lnpn = WEIGHTS .* lnpn; 
f = -1 * sum(lnpn);

	
			
