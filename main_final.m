% This program compute the MLE estimator for binary mixed logit model with 
% panel data
clear all;


% Define global variables
global REPTITION NOBS NG XMAT YVEC TIMES NVAR NFC IDFC NNC IDNC...
	NUC IDUC NLC IDLC NERROR UERROR LERROR REVEAL STATED CALPOLY... 
	BROOK TRANS CVAR SERROR;  

% REPTITION is the number of random draws   
REPTITION = 1500;

% The way of draw
HALTON = 1;

% NOBS is the number of observations    
NOBS=1155;

% NG is the number of individuals
NG=548;

% XMAT is the matrix containing all explanatory variables		
%%%%%%%%%%%%%%%%%%%%%%%%%
load cal_broo_sp_test1.txt;
data=cal_broo_sp_test1;
clear cal_broo_sp_test1;
id = data(:,1);
route = data(:,2);
hsize = data(:,3);
dist = data(:,4);
female = data(:,5);
median = data(:,6);
dmp80 = data(:,7);
cost = data(:,8);
dlate = data(:,9);
hc = data(:,10);
mc = data(:,11);
dage = data(:,12);
calpoly = data(:,13);
brook = data(:,14);
sst = data(:,15);
llt = data(:,16);
rp = data(:,17);
sp = data(:,18);
date = data(:,19);
group = data(:,20);
one = ones(NOBS,1);
clear data;

d2=dist.^2;
d3=dist.^3;
dt=dist.*median.*rp;
dt2=d2.*median.*rp;
dt3=d3.*median.*rp;
sst=sst.*sp;
llt=llt.*sp;

cmedian=calpoly.*median;
bmedian=brook.*median;
rpmedian=rp.*median;
spmedian=sp.*median;

cdmp80=calpoly.*dmp80;
bdmp80=brook.*dmp80;
rpdmp80=rp.*dmp80;
spdmp80=sp.*dmp80;

rpcost=rp.*cost;
spcost=sp.*cost;

rphc=rp.*hc;
sphc=sp.*hc;

rpmc=rp.*mc;
spmc=sp.*mc;

female_rp = female .* rp;
female_sp = female .* sp;

dlate_rp = dlate .* rp;
dlate_sp = dlate .* sp;

dage_rp = dage .* rp;
dage_sp = dage .* sp;

hsize_rp = hsize .* rp;
hsize_sp = hsize .* sp;

%%%%%%%%%%%%%%%%%%%%%%%%%
REVEAL=rp;
STATED=sp;
CALPOLY=calpoly;
BROOK=brook;

XMAT=[sp,median,dmp80,one,brook,calpoly,rpcost,rphc,rpmc,dt,dt2,dt3,...
rpdmp80,spcost,sphc,spmc,sst,llt,spdmp80,female,dlate,dage,hsize]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NOBS;
	if date(i)~=1
		group(i)=0;
	end;
end;
clear i;

% YVEC is the vector of dependent variable
YVEC = route;
%clear id route median cost dmp75 female nincome dlate one sp;

% TIMES is a NG X 1 vector, which indicates the number of choice occassion
% each individual faces
if NG~=NOBS
	TIMES=nonzeros(group);
else
	TIMES=ones(NOBS,1);
end;

if sum(TIMES)~=NOBS
	display wrong;
end;

% NVAR is the number of explanatory variables
NVAR=23;

% NFC is the number of varibles with fixed coefficients
NFC=19;

% IDFC indicates which variables have fixed coefficients
IDFC=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];

% NNC is the number of variables with normal distributed coefficients
NNC=4;

% IDNC indicates which variables have normal coeffiicents
IDNC=[1,2,3,4];    

% Generate Normal Random Numbers
if NNC>0 & HALTON == 0
	randn('seed', 1234)
  	ndraw=randn(NG,NNC*REPTITION);
	if NOBS==NG
		NERROR=ndraw;
	else
	     	NERROR=zeros(NOBS,NNC*REPTITION);
        	ini=0;
        	for j=1:NG;
        		gstart=ini+1;
        		gend=ini+TIMES(j);
        		temp0=ones(TIMES(j),1);
        		temp1=kron(ndraw(j,:),temp0);
        		NERROR(gstart:gend,:)=temp1;
        		ini=ini+TIMES(j);
        		clear gstart gend temp0 temp1;
    		end;
	end;
else
	NERROR=[0];
end;
clear j ndraw;

% Generate Normal Random Numbers Using Halton Sequences
% Pime numbers to generate halton sequences
PRIME = [11;89;19;23;29;31;37;41;43;47;53;59;61;67;71;73;79;83;89;97;101];
if NNC>0 & HALTON == 1
	NERROR = zeros(sum(TIMES), NNC*REPTITION);
	temp = zeros(NG, NNC*REPTITION); 
	for i = 1:NNC;
		gstart = (i-1)*REPTITION + 1;
		gend = (i-1) * REPTITION + REPTITION;
      		seed = PRIME(i);
       		hdraw = halton( NG*REPTITION, 0, seed, 1);
       		hdraw = transpose(hdraw);
       		for j = 1:NG;
       			temp(j,gstart:gend) =...
            		hdraw(1,((j-1)*REPTITION + 1) : ((j-1)*REPTITION+REPTITION) );
       		end;
		clear hdraw;
	end;
	clear gstart gend;
       		
	ini=0;
      	for i = 1:NG;
       		gstart = ini+1;
       		gend = ini+TIMES(i);
       		NERROR(gstart:gend,:) = kron( temp(i,:), ones(TIMES(i),1) );
       		ini=ini+TIMES(i);
       	end;
       	clear temp;
end;


% Set initial values for parameters
%b0 =[-1;-1;-0.6;0.3;0.1;-0.1;0.0;0.0;-0.2;-0.5;0.0;0.0;-0.1;;-0.1;-2;...
%0.4;0.1;0.0;-0.1;0;0;0;0;-1;1;0.0;0.0;0.0;1.0]; 
%scale1=[1];
%scale2=[1];
%b=[b0;scale1;scale2];
load betanew_r3000.txt;
b = betanew_r3000;
f = simullikelihood(b)
%STOP
%Set options
options=optimset('TolFun',1*10^-2,'Display','iter','MaxFunEvals',300000000); 

% Find maximum of simulated log-likelihood function 
[beta,fval,exitflag,output,grad,hessian]=fminunc('simullikelihood',... 
b, options);
stop
% The following part is to calculate robust std. error numerically
deriv = zeros( NG, length(beta) );
step_length = 1*10^-6;
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

