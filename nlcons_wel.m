function [c, ceq] = nlcons_wel(toll)

global KE KF FIXEDE FIXEDF CONSTRUCTION DISCOUNT NPOP BPOWER;

%tolle = toll(1);
%tollf = toll(2);

junk = equilibrium(toll, FIXEDF, 3);

veupdate = sum(junk(:,1) + junk(:,3)/2 + junk(:,5)/3);
vfupdate = sum(junk(:,2) + junk(:,4)/2 + junk(:,6)/3);


te = supply(KE, veupdate);
tf = supply(KF, vfupdate);
re = (te - 9.23) * 0.3785;
rf = (tf - 9.23) * 0.3785;

%profit = (tollf * vfupdate + tolle * veupdate) * 0.65 * 250;
%profit = toll * veupdate * 0.65 * 250;
%f = (profit / DISCOUNT) -  (CONSTRUCTION);
%c(1) = -1 * f;

%% Calculate social welfare under equilibrium %%
[wel] = welfare(toll, FIXEDF, te, tf, re, rf, 3);
f= (wel/NPOP) - (12.12);
c = -1 * f;
ceq = [];
