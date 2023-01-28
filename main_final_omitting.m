% This program compute the MLE estimator for binary mixed logit model with 
% panel data
clear all;


% Define global variables
global REPTITION NOBS NG XMAT YVEC TIMES NVAR NFC IDFC NNC IDNC...
	NUC IDUC NLC IDLC NERROR UERROR LERROR REVEAL STATED CALPOLY... 
	BROOK TRANS CVAR SERROR SPONLY RSP;  

% REPTITION is the number of random draws   
REPTITION=1000;

% NOBS is the number of observations    
NOBS=1155;

% NG is the number of individuals
NG=548;

% XMAT is the matrix containing all explanatory variables		
%%%%%%%%%%%%%%%%%%%%%%%%%
load econometrica.txt;
data = econometrica;
clear econometrica;
id = data(:,1);
route = data(:,2);
hsize = data(:,3);
dist = data(:,4);
ptime = data(:,5);
dage = data(:,6);
female = data(:,7);
median = data(:,8);
dmp80 = data(:,9);
cost = data(:,10);
dlate = data(:,11);
commute = data(:,12);
hc = data(:,13);
mc = data(:,14);
ll = data(:,15);
ss = data(:,16);
calpoly = data(:,17);
brook = data(:,18);
sp = data(:,19);
rsp = data(:,20);
date = data(:,21);
group = data(:,22);
one = ones(NOBS,1);
clear data;
for i = 1: NOBS;
	if id(i) == 1235782 
		dage(i) = 1;
	end;
end;
	

sponly = zeros(NOBS,1);
peak = zeros(NOBS,1);
rp = zeros(NOBS,1);
for i = 1:NOBS;
	if rsp(i) == 0
		sponly(i) = 1;
	end;
	if ptime(i)>=700 & ptime(i)<=800
		peak(i) = 1;
	end;
	if calpoly(i) == 1 | brook(i) == 1
		rp(i) = 1;
	end;
end;
sponly = sponly .* sp;
rsboth = brook + rsp;

d2 = dist.^2;
d3 = dist.^3;
dt = dist.*median.*rp;
dt2 = d2.*median.*rp;
dt3 = d3.*median.*rp;
st = ss.*median.*sp;
lt = ll.*median.*sp;

rpdmp80 = rp.*dmp80;
spdmp80 = sp.*dmp80;

rpcost = rp.*cost;
spcost = sp.*cost;

rphc = rp.*hc;
sphc = sp.*hc;

rpmc = rp.*mc;
spmc = sp.*mc;

%%%%%%%%%%%%%%%%%%%%%%%%%
REVEAL = rp;
STATED = sp;
CALPOLY = calpoly;
BROOK = brook;
SPONLY = sponly;
RSP = rsp;

XMAT=[sp,cost,median,dmp80,one,brook,calpoly,rpcost,rphc,rpmc,dt,...
dt2,dt3,rpdmp80,spcost,sphc,spmc,st,lt,spdmp80,dage,hsize]; 
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
NVAR=22;

% NFC is the number of varibles with fixed coefficients
NFC=17;

% IDFC indicates which variables have fixed coefficients
IDFC=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];

% NNC is the number of variables with normal distributed coefficients
NNC=5;

% IDNC indicates which variables have normal coeffiicents
IDNC=[1,2,3,4,5];    

% NUC is the number of variables with uniform distributed coefficients
NUC=0;

% IDUC indicates which variables have uniform coefficients
IDUC=[0];
 
% NLC is the number of variables with log-normal distributed coefficients
NLC=0;

% IDLC indicates which variables have log-normal coefficients
IDLC=[0];

% CVAR is the # of  variables which are constrained to be zero  
CVAR=0;

% IDCVAR indicates which variables are constrained to be zero
IDCVAR=[0];

% Generate Normal Random Numbers
if NNC>0
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

% Generate Uniform Random Numbers
if NUC>0
        udraw=((rand(NG,NUC*REPTITION)*2)-1)*3;
        UERROR=zeros(NOBS,NUC*REPTITION);
        ini=0;
        for j=1:NG;
        	gstart=ini+1;
          	gend=ini+TIMES(j);
                temp0=ones(TIMES(j),1);
                temp1=kron(udraw(j,:),temp0);
                UERROR(gstart:gend,:)=temp1;
                ini=ini+TIMES(j);
                clear gstart gend temp0 temp1;
         end;
         clear j;
else
	UERROR=[0];
end;
clear udraw;                   

% genertae normal error for log-normal distributed coefficients      
if NLC>0
        ldraw=randn(NG,NLC*REPTITION);
        LERROR=zeros(NOBS,NLC*REPTITION);
        ini=0;
        for j=1:NG;
                gstart=ini+1;
                gend=ini+TIMES(j);
                temp0=ones(TIMES(j),1);
                temp1=kron(ldraw(j,:),temp0);
                LERROR(gstart:gend,:)=temp1;
                ini=ini+TIMES(j);
                clear gstart gend temp0 temp1;
        end;
else
        LERROR=[0];
end;
clear j ldraw;                                                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SERROR=randn(NOBS,REPTITION);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for parameters
b0 = [0;
	-2;
	-2;
	1;
	0.6;
	-0.4;
	0.1;
	-0.004;
	-0.5;
	-1;
	0;
	0;
	-0.2;
	-0.2;
	-5;
	1;
	-0.4;
	-1;
	0;
	0;
	0;
    1;
    1];
scale1=[1];
scale2=[1];
b=[b0;scale1;scale2];

% Set options
options=optimset('TolFun',1*10^-3,'Display','iter','MaxFunEvals',300000000); 

% Find maximum of simulated log-likelihood function 
[beta,fval,exitflag,output,grad,hessian]=fminunc('simullikelihood',... 
b, options);

