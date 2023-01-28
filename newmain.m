% This program compute the MLE estimator for mixed logit model combining 
% RP and SP data and with a unbalance panel data 
% Jia Yan, 11/07/2003
%clear all;

% Define global variables
global REPTITION NUSERS XFIXED XRANDOM YVEC NFC NNC SCALE...
       CALPOLY SP NALT WEIGHTS HALTON NERROR SCALE_CAL SCALE_BR SCALE_SP...
       BROOK TIMES CHOICE RP;  


% Do you want to use halton random draws? 1 if yes, 0 then use quasi-random generator.
HALTON = 1;

% Choosing the way of scaling variance of different data set
% 1 if adding adidtional random term to data set with smaller variance
% 0 if multiplying the equation with scale parameter
% sometimes, it is easier to get convergent results using 1   
SCALE = 0;

% REPTITION is the number of random draws used to approximate integration   
REPTITION = 4000;

% NUSERS is the number of persons
NUSERS = 538;

% NALT is the number of RP choice alternatives
NALT = 9;

% NFC is the number of varibles with fixed coefficients
NFC = 26;

% NNC is the number of variables with normal distributed coefficients
NNC = 5;

%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XMAT is the matrix containing all explanatory variables.
% YVEC is the vector of dependent variable.
% Combining different data sets is feasible. You need specify the different
% data sources, like CALPOLY here.
% Data can be generated from statistics package such as STATA.
%%%%%%%%%%%%%%%%%%%%%%%%%
load rpspdata.txt;
rpspdata = rpspdata;
CHOICE = rpspdata(:,2);
YVEC = rpspdata(:,3);
XFIXED = rpspdata(:,4:29);
XRANDOM = [rpspdata(:,30:34)];
BROOK = rpspdata(:,35);
CALPOLY = rpspdata(:,36);
SP = XFIXED(:,26);
RP = BROOK + CALPOLY;
clear rpspdata;

% WEIGHTS is used for choice-based sample 
load rpspweights.txt;
WEIGHTS = rpspweights;

% TIMES is NUSERS by 1 vector recording number of observations for each person 
load rpsptimes.txt;
TIMES = rpsptimes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRIME = [2;3;5;7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;79;83;89;97;101];

% Generate Normal Random Numbers Using Halton Sequences
if NNC>0 & HALTON == 1
	NERROR = zeros(sum(TIMES), NNC*REPTITION); 
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
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(TIMES(i),1) );
       		ini=ini+TIMES(i);
       	end;
       	clear temp;
   
end;

if NNC > 0 & HALTON == 0
	%randn('state', 13457);
    	NERROR = zeros( NUSERS, NNC*REPTITION);
    	temp = randn(NUSERS, NNC*REPTITION);
    	ini=0;
      	for i = 1:NUSERS;
       		gstart = ini+1;
       		gend = ini+TIMES(i);
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(TIMES(i),1) );
       		ini=ini+TIMES(i);
       	end;
       	clear temp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if HALTON == 0
	SCALE_SP = randn(sum(TIMES), REPTITION);
	SCALE_CAL = randn(sum(TIMES), REPTITION);
	SCALE_BR = randn(sum(TIMES),REPTITION);
	for i = 1: sum(TIMES);
		if CHOICE(i)==3 | CHOICE(i) == 6 | CHOICE(i) == 9
       			SCALE_BR(i,:) = 0;    
       		end;
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load beta_revised_nhet.txt;
b = beta_revised_nhet;

%b0 =[b(1:27); 0.1; b(28:31); 1; b(32:33)];
b0 = b;

%% randomization to check the robustness %%
b0(15:30) = b0(15:30) + randn(16,1);

% Set options
options=optimset('TolFun',1*10^-1,'Display','iter','MaxFunEvals',300000000); 

% Find maximum of simulated log-likelihood function 
[beta,fval,exitflag,output,grad,hessian]=fminunc('newlikelihood',... 
b0, options);

% The following part is to calculate robust std. error numerically
deriv = zeros( NUSERS, length(beta) );
step_length = 1*10^-4;
for i = 1:length(beta);
	nbeta1 = beta;
	nbeta2 = beta;
	nbeta1(i,1) = nbeta1(i,1) + step_length/2;
	nbeta2(i,1) = nbeta2(i,1) - step_length/2;
	deriv(:,i) = (func(nbeta1) - func(nbeta2)) ./ (step_length);  
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
