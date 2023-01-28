
nperiods = 100;

w = 1/0.999995;
c = 0.7;
pv = 0.99;  

yini = 0.79;
res = zeros(nperiods, 2);
for i = 1: nperiods;
	i
	yupdate = (w / yini) * (pv - c) + pv - w;
	if yupdate > 1
		yupdate = 1;
	end;
	res(i,1) = i;
	res(i,2) = yupdate;
	yini = yupdate;
end;

	