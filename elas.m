% This program compute the MLE estimator for mixed logit model combining 
% RP and SP data and with a unbalance panel data 
% Jia Yan, 11/07/2003
clear all;

% Define global variables
global REPTITION NUSERS XFIXED XRANDOM YVEC NFC NNC SCALE...
       CALPOLY SP NALT WEIGHTS HALTON NERROR SCALE_CAL SCALE_BR SCALE_SP...
       BROOK TIMES CHOICE RP;  


% Do you want to use halton random draws? 1 if yes, 0 then use quasi-random generator.
HALTON = 1;

% REPTITION is the number of random draws used to approximate integration   
REPTITION = 3000;

% NUSERS is the number of persons
NUSERS = 79;

% NALT is the number of RP choice alternatives
NALT = 9;

% NFC is the number of varibles with fixed coefficients
NFC = 26;

% NNC is the number of variables with normal distributed coefficients
NNC = 7;

PRIME = [2;3;5;7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;79;83;89;97;101];

% Generate Normal Random Numbers Using Halton Sequences
if NNC>0 & HALTON == 1
	NERROR = zeros(NALT*NUSERS, NNC*REPTITION); 
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
       		gend = ini+NALT;
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(NALT,1) );
       		ini=ini+NALT;
       	end;
       	clear temp;
   
end;

%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load beta_revised_het.txt;
beta = beta_revised_het * beta_revised_het(34,1);
beta(23:25) = beta(23:25) / beta_revised_het(34,1);
beta(28) = beta(28) / beta_revised_het(34,1);
beta(33) = beta(33) / beta_revised_het(34,1);
for i = 1:2;
	if i == 1
		load brdata0.txt;
		rpspdata = brdata0;
	else
		load brdata1.txt;
		rpspdata = brdata1;
	end;
	ID = rpspdata(:,1);
	CHOICE = rpspdata(:,2);
	YVEC = rpspdata(:,3);
	XFIXED = rpspdata(:,4:29);
	XRANDOM = [rpspdata(:,30:36)];
	BROOK = rpspdata(:,37);
	CALPOLY = rpspdata(:,38);
	SP = XFIXED(:,26);
	RP = BROOK + CALPOLY;
	clear rpspdata;

	temp = marginal(beta);

	if i == 1
		res = [ID, CHOICE, temp];
	else
		res = [res, temp];
	end;
end;

save fjunk.txt res -ASCII;


