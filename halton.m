function [phi] = halton(m, rw, s, order)
% This function generates standard normal random numbers based on halton sequence.
% The halton sequence here are randomized and shuffled to avoid correlation
% in high dimension integration. 
% The two parameters are, m: length of generated numbers;
% rw: integer to randomized halton sequences. The logic here is first
%    to draw an integer from 0 to k with equal probability, and label this 
%    integer as discard; the generated halton sequence is resulted by discarding
%    the first "discard" numbers;  s : prime number used; 
% 
% Jia Yan, 11/08/2003

%% First step: draw an integer between zero and rw %%
if rw > 0 
	uniform = rand(1,1);
	for i = 1:(k+1);
		low = (i-1)/(k+1);
		upper = i/(k+1);
		if uniform > low & uniform <= upper
			discard = i - 1;
			break;
		end; 
	end;
else 
	discard = 0;
end;

%% After above randomnization, each sequence will drop the first "discard" numbers %% 		     	
n = m + discard;                   
k=fix(log(n+1)/log(s));      
phi=[0];
i=1;
while i<=k,
	x=phi;
	j=1;
	while j<s,
		y=phi+(j/s^i);
		x=[x;y];
		j=j+1;
	end;
	phi=x;
	i=i+1;
end;

x=phi;
j=1;
while j<s & length(x) < n+1,
	y=phi+(j/s^i);
	x=[x;y];
	j=j+1;
end;

phi=x( (2+discard):(n+1) );


% Randomized by adding each element with a uniform r.v. %%
if rw == 0 
	urand = rand(1,1);
	for i = 1:length(phi);
		phi(i) = phi(i) + urand;
		if phi(i) >= 1
			phi(i) = phi(i) - 1;
		end;
	end;
	
	temp = nonzeros(phi);
	if length(temp) ~= length(phi)
		display('error')
	end;
end;

%% Now, we have obtained randomized halton sequences. To avoid the correlation in %%
%% high dimension problem, we shuffle the sequence by reordering it based on the rank of %%
%% elements in a uniform random vector "uniform". %%
if order == 1
	uniform = rand(m,1);
	[uniform,ind] = sort(uniform);
	phi = phi(ind);
end;

phi = norminv(phi,0,1);
%for i = 1:length(phi);
%	phi(i) = norminv(phi(i), 0, 1);
%end;
   


