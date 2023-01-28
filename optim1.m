
function [f] = optim1(toll)
% This function is for one-route pricing case -- to get travel times, optimal toll,
% as well as welfare under differnt toll levels on express lanes. 

global KE KF FIXEDE FIXEDF;  

[junk] = equilibrium(toll, FIXEDF, 3);

veini = junk(:,1) + junk(:,3)/2 + junk(:,5)/3;
vfini = junk(:,2) + junk(:,4)/2 + junk(:,6)/3;	
veupdate = sum(veini);
vfupdate = sum(vfini);

te = supply(KE, veupdate);
tf = supply(KF, vfupdate);
re = (te - 9.23) * 0.3785;
rf = (tf - 9.23) * 0.3785;

%% Calculate social welfare under equilibrium %%
[wel] = welfare(toll, FIXEDF, te, tf, re, rf, 3);
%f =  veupdate * toll - 0.0008 * 10 * veupdate;
%f = -1 * f;

f = -1 * wel;

