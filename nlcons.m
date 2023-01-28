
function [c, ceq] = nlcons(toll)

global KE KF FIXEDE FIXEDF CONSTRUCTION DISCOUNT;

%tolle = toll(1);
%tollf = toll(2);

junk = equilibrium(toll, FIXEDF, 3);

veupdate = sum(junk(:,1) + junk(:,3)/2 + junk(:,5)/3);
vfupdate = sum(junk(:,2) + junk(:,4)/2 + junk(:,6)/3);


%profit = (tollf * vfupdate + tolle * veupdate) * 0.65 * 250;
profit = toll * veupdate * 0.65 * 250;
f = (profit / DISCOUNT) - (CONSTRUCTION)*(2/3);

c = -1 * f;
%c = [];
ceq = [];

