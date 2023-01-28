function [junk, counter] = equilibrium(tolle, tollf, s)

%%%%% Search Equilibrium %%%%%%
% Define global variables
global NALT KE KF STEPLENGTH; 

%% Set initial traffic conditions %%
ini = 100 * ones(NALT,1); 

%% Update %%
bench = supply(KE,0);
updated = zeros(NALT,1);
cr = 1;	
counter = 0;
while cr == 1,
	counter = counter + 1;
		
	if counter > 1000 
		display('no convergence');
		break;
	end;
	
	veini = ini(1) + ini(3)/2 + ini(5)/3;
	vfini = ini(2) + ini(4)/2 + ini(6)/3;
	teupdate = supply(KE, veini);
	tfupdate = supply(KF, vfini);
	reupdate = (teupdate - bench) * 0.3785;
	rfupdate = (tfupdate - bench) * 0.3785;
	
	[h, junk] = demand(tolle, tollf, teupdate, tfupdate, reupdate, rfupdate, s);
	cr = 0;
	for i = 1:NALT;
		diff = sum(junk(:,i)) - ini(i);
		updated(i,1) = diff; 
		ini(i) = ini(i) + STEPLENGTH * diff;
	
		if abs(diff) > 5 & cr == 0
			cr = 1;	
		end;
	end;
end;
