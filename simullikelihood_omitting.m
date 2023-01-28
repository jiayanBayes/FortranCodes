function [sl] = simullikelihood(beta)
% This function calculates simulated likelihood function for mixed logit
% model
% Declare global variables
global REPTITION NOBS NG XMAT YVEC TIMES NVAR NFC IDFC NNC IDNC...
	NUC IDUC NLC IDLC NERROR UERROR LERROR REVEAL STATED TRANS...
	CVAR IDCVAR SERROR CALPOLY BROOK RSP SPONLY; 

% Calculate the conditional probability for each individual and each
% time periods

pnt=zeros(NOBS,REPTITION);
for i=1:REPTITION;
	x=zeros(NOBS,(NFC+2*NNC+2*NUC));
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
		iter1=1;
		%iter2=(i-1)*NNC+1;
		for j=1:NNC;
			nstart=(j-1)*REPTITION+i;
			x(:,NFC+iter1)=XMAT(:,IDNC(j));
			x(:,NFC+iter1+1)=NERROR(:,nstart).*XMAT(:,IDNC(j));
			%x(:,NFC+iter1+1)=NERROR(:,iter2).*XMAT(:,IDNC(j));
			%iter2=iter2+1;
			iter1=iter1+2;
		end;	
		clear j;
	else
		x=x;
	end;
	clear iter1 iter2;
	
	% Add uniform random terms into XMAT
 	if NUC>0
		iter1=1;
                iter2=(i-1)*NUC+1;
                for j=1:NUC;
			x(:,NFC+2*NNC+iter1)=XMAT(:,IDUC(j));
                        x(:,NFC+2*NNC+iter1+1)=UERROR(:,iter2).*...
				XMAT(:,IDUC(j));
                        iter2=iter2+1;
			iter1=iter1+2;
                end;
                clear j;
	else
		x=x;
        end;
        clear iter1 iter2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	x1 = x(:,1:19);
	x2 = x(:,21);
	x3 = x(:,23);
	x4 = x(:,25);
	x5 = x(:,27);
	nx = [x1,x2,x3,x4,x5];      
 	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	index=nx(:,1:21) * beta(1:21);
	for j=1:NOBS;
		if REVEAL(j) == 1
			index(j) = index(j) + beta(9) * beta(22) * ...
				nx(j,22) + nx(j,23);
		end;
		if STATED(j) == 1
			index(j) = index(j) + beta(15) * beta(22) *... 
				nx(j,22) + beta(23) * nx(j,23);
		end;
    end;
	clear j;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% index here includes all the fixed, normal and uniform variables
	% rescaling
	for j = 1:NOBS;
       	if BROOK(j) == 1 
           	index(j) = index(j);
		end;
		if CALPOLY(j) == 1
			index(j) = beta(24) * index(j);
		end;		
		if STATED(j) == 1
		    index(j) = beta(25) * index(j);
       	end;
    end;
	clear j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	index = exp(-1*index);
	for j = 1:NOBS;
		if YVEC(j)==1
			pnt(j,i)=1/(1+index(j));
		else
			pnt(j,i)=index(j)/(1+index(j));
		end;
	end;
	clear j
	clear index x;
end;
clear i;

% Calculate the joint probability of a sequence of choices made by same
% individuals
pn = zeros(NG,REPTITION);
spn = zeros(NG,1);
ini = 0;
for i = 1:NG;
	gstart = ini+1;
	gend = ini+TIMES(i);
	temp = pnt(gstart:gend,:);
	for j=1:REPTITION;
		pn(i,j) = prod(temp(:,j));
	end;
	spn(i)=(sum(pn(i,:)))/REPTITION;
	clear j gstart gend temp;
	ini=ini+TIMES(i);
end;
clear i;


lnspn=log(spn);
sl=-1*sum(lnspn);

	
			
