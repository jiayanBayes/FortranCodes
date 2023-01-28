% This program compute the MLE estimator for mixed logit model combining 
% RP and SP data and with a unbalance panel data 
% Jia Yan, 05/26/2003, ChoiceStream Inc.
clear all;

% Define global variables
global REPTITION NUSERS XMAT YVEC NFC IDFC NNC IDNC SCALE...
	CALPOLY	SP NALT WEIGHTS HALTON HRANDN HRANDN_SCALE_CAL...
	HRANDN_SCALE_SP HRANDN_SCALE_BR NERROR BROOK TIMES CHOICE RP;  


% Do you want to use halton random draws? 1 if yes, 0 then use quasi-random generator.
HALTON = 0;

% Choosing the way of scaling variance of different data set
% 1 if adding adidtional random term to data set with smaller variance
% 0 if multiplying the equation with scale parameter
% sometimes, it is easier to get convergent results using 1   
SCALE = 1;

% REPTITION is the number of random draws used to approximate integration   
REPTITION = 800;

% NUSERS is the number of persons
NUSERS = 538;

% NALT is the number of RP choice alternatives
NALT = 9;

% XMAT is the matrix containing all explanatory variables.
% YVEC is the vector of dependent variable.
% Combining different data sets is feasible. You need specify the different
% data sources, like CALPOLY here.
% Data can be generated from statistics package such as STATA.
%%%%%%%%%%%%%%%%%%%%%%%%%
load rpspdata.txt;
CHOICE = rpspdata(:,2);
YVEC = rpspdata(:,3);
XMAT = rpspdata(:,4:46);
BROOK = rpspdata(:,47);
CALPOLY = rpspdata(:,48);
SP = rpspdata(:,4);
RP = BROOK + CALPOLY;
clear rpspdata;

% WEIGHTS is used for choice-based sample 
load rpspweights.txt;
WEIGHTS = rpspweights;

% TIMES is NUSERS by 1 vector recording number of observations for each person 
load rpsptimes.txt;
TIMES = rpsptimes;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% NFC is the number of varibles with fixed coefficients
NFC = 34;

% IDFC indicates which variables have fixed coefficients
IDFC = [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,...
	29,30,31,32,33,34,35,36,37,38,39,40,41,42,43];

% NNC is the number of variables with normal distributed coefficients
NNC = 9;

% IDNC indicates which variables have normal coeffiicents
IDNC = [1,2,3,4,5,6,7,8,9];    

% Pime numbers to generate halton sequences
PRIME = [17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;79;83;89;97;101;2;3;5;7;11];

% Generate Normal Random Numbers Using Halton Sequences
if NNC>0 & HALTON == 1
	HRANDN = zeros(sum(TIMES), NNC*REPTITION);
	temp = zeros(NUSERS, NNC*REPTITION); 
	for i = 1:NNC;
		gstart = (i-1)*REPTITION + 1;
		gend = (i-1) * REPTITION + REPTITION;
      		seed = PRIME(i);
       		hdraw = halton( NUSERS*REPTITION, 0, seed, 1);
       		hdraw = transpose(hdraw);
       		for j = 1:NUSERS;
       			temp(j,gstart:gend) =...
            		hdraw( 1, ((j-1)*REPTITION + 1) : ((j-1)*REPTITION+REPTITION) );
       		end;
		clear hdraw;
	end;
	clear gstart gend;
       		
	ini=0;
      	for i = 1:NUSERS;
       		gstart = ini+1;
       		gend = ini+TIMES(i);
       		HRANDN(gstart:gend,:) = kron( temp(i,:), ones(TIMES(i),1) );
       		ini=ini+TIMES(i);
       	end;
       	clear temp;
end;

if NNC > 0 & HALTON == 0
	NERROR = zeros(sum(TIMES), NNC*REPTITION);
	randn('state', 13457);
    	ndraw = randn( NUSERS, NNC*REPTITION);
    	ini=0;
        for i = 1:NUSERS;
        	gstart = ini+1;
        	gend = ini+TIMES(i);
        	NERROR(gstart:gend,:) = kron( ndraw(i,:), ones(TIMES(i),1) );
        	ini=ini+TIMES(i);
        end;
    	clear ndraw;
end;


% HRANDN_SCALE if halton random draws to resacle variance of different data set 
if SCALE == 1 & HALTON == 1
	HRANDN_SCALE_CAL = zeros(sum(TIMES), REPTITION);
	HRANDN_SCALE_SP = zeros(sum(TIMES), REPTITION);
	HRANDN_SCALE_BR = zeros(sum(TIMES), REPTITION);
	hdraw_scale_cal = halton( sum(TIMES)*REPTITION, 400, PRIME(NNC+1));
	hdraw_scale_sp = halton( sum(TIMES)*REPTITION, 400, PRIME(NNC+2));
	hdraw_scale_br = halton( sum(TIMES)*REPTITION, 400, PRIME(NNC+3));
	for j = 1:NUSERS;
   		HRANDN_SCALE_CAL(j,:) = ...
   		hdraw_scale_cal( ((j-1)*REPTITION + 1):((j-1)*REPTITION+REPTITION));
   		HRANDN_SCALE_SP(j,:) = ...
   	        hdraw_scale_sp( ((j-1)*REPTITION + 1):((j-1)*REPTITION+REPTITION));
   	        HRANDN_SCALE_BR(j,:) = ...
   	        hdraw_scale_br( ((j-1)*REPTITION + 1):((j-1)*REPTITION+REPTITION));
	end;
	clear hdraw_scale_cal hdraw_scale_sp hdraw_scale_br;
    
    for j = 1: sum(TIMES);
        if CHOICE(j)==3 | CHOICE(j) == 6 | CHOICE(j) == 9
            HRANDN_SCALE_BR(j,:) = 0;    
        end;
    end;
end;

if SCALE == 1 & HALTON == 0 
	HRANDN_SCALE_SP = randn(sum(TIMES),REPTITION);
	HRANDN_SCALE_CAL = randn(sum(TIMES), REPTITION);
	HRANDN_SCALE_BR = randn(sum(TIMES), REPTITION);
	for j = 1: sum(TIMES);
        	if CHOICE(j)==3 | CHOICE(j) == 6 | CHOICE(j) == 9
            		HRANDN_SCALE_BR(j,:) = 0;    
        	end;
    	end;
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for parameters
b_fixed =[-1;
	0.5;
	0.3;
    	-0.2;
    	0.04;   
    	-0.001;
    	-0.5;
    	-0.4;
    	-0.4;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
    	0;
        0;
        1;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0];
        
b_ran = [0;
	0;
	0;
	0;
	0;
	0;
	0];	
s=[1;1];
b=[b_fixed;b_ran;s];

% Set options
options=optimset('TolFun',1*10^-3,'Display','iter','MaxFunEvals',300000000); 

% Find maximum of simulated log-likelihood function 
[beta,fval,exitflag,output,grad,hessian]=fminunc('mixed_logit_rpsp',... 
b, options);

% The following part is to calculate robust std. error numerically
deriv = zeros( NUSERS, length(beta) );
step_length = 1*10^-6;
for i = 1:length(beta);
	nbeta1 = beta;
	nbeta2 = beta;
	nbeta1(i,1) = nbeta1(i,1) + step_length;
	nbeta2(i,1) = nbeta2(i,1) - step_length;
	deriv(:,i) = (likelihood_func(nbeta1) - likelihood_func(nbeta2)) ./ (2*step_length);  
end;	

BHHH = deriv' * deriv;
H = inv(hessian);
robust = H * BHHH;
robust = robust * H;

std_h = diag(H);
std_h = sqrt(std_h);
std_bhhh = inv(BHHH);
std_bhhh = diag(std_bhhh);
std_bhhh = sqrt(std_bhhh);
std_rob = diag(robust);
std_rob = sqrt(std_rob);
  
RESULTS = [beta, std_h, std_bhhh, std_rob]
