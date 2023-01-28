
function [f] = optim2(toll)
% This function is for two-route pricing case -- to get travel times, optimal toll,
% as well as welfare under differnt toll levels on express lanes. 

global NPOP CONS RHO PARAMETER KE KF DISCOUNT CONSTRUCTION WEIGHT;  

tolle = toll(1,1);
tollf = toll(2,1);

[junk] = equilibrium(tolle, tollf, 3);

veini = junk(:,1) + junk(:,3)/2 + junk(:,5)/3;
vfini = junk(:,2) + junk(:,4)/2 + junk(:,6)/3;	
veupdate = sum(veini);
vfupdate = sum(vfini);

te = supply(KE, veupdate);
tf = supply(KF, vfupdate);
re = (te - 9.23) * 0.3785;
rf = (tf - 9.23) * 0.3785;

%% Calculate social welfare under equilibrium %%
[wel] = welfare(tolle, tollf, te, tf, re, rf, 3);

profit = tolle * veupdate + tollf * vfupdate - 0.0008 * 10 * (veupdate + vfupdate);
%f = -1 * profit;
f = wel;
%f = wel + profit;
%f = -1 * f;

%profit = (tollf * vfupdate + tolle * veupdate) * 0.65 * 250;
%profit = (profit / DISCOUNT) -  CONSTRUCTION;

%f(1) =  wel;
%f(2) = profit;

	 