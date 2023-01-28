
function [f] = optim1_sequential(toll)
% This function is for one-route pricing case -- to get travel times, optimal toll,
% as well as welfare under differnt toll levels on express lanes. 

global KE KF DISCOUNT CONSTRUCTION;  

load Response_Bargaining_Profit_Part1.txt;
load Response_Bargaining_Profit_Part2.txt;

if toll > 3.52 & toll <= 3.58
	toll = 3.51;
end;

if toll > 3.58 & toll <= 3.64
	toll = 3.64;
end;


if toll <= 3.58
	x = Response_Bargaining_Profit_Part1(:,2);
	y = Response_Bargaining_Profit_Part1(:,1);
else
	x = Response_Bargaining_Profit_Part2(:,2);
	y = Response_Bargaining_Profit_Part2(:,1);
end;

tolle = interp1(x, y, toll)

[junk] = equilibrium(tolle, toll, 3);

veini = junk(:,1) + junk(:,3)/2 + junk(:,5)/3;
vfini = junk(:,2) + junk(:,4)/2 + junk(:,6)/3;	
veupdate = sum(veini);
vfupdate = sum(vfini);

te = supply(KE, veupdate);
tf = supply(KF, vfupdate);
re = (te - 9.23) * 0.3785;
rf = (tf - 9.23) * 0.3785;

%% Calculate social welfare under equilibrium %%
[wel] = welfare(tolle, toll, te, tf, re, rf, 3);
%f =  vfupdate * toll - 0.0008 * 10 * vfupdate;
%f = -1 * f;

f(1) = wel;

profit = toll * vfupdate * 0.65 * 250;
f(2) = (profit / DISCOUNT) - 0.5 * CONSTRUCTION;

