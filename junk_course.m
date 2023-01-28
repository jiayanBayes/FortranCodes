
obs = 1000; 
errors = randn(obs, 1);

y1 = zeros(obs,1);
y2 = zeros(obs,1);
y3 = zeros(obs,1);

y4 = zeros(obs,1);
y5 = zeros(obs,1);
y6 = zeros(obs,1);
for i = 1:obs;
	if i == 1
		y1(i,1) = randn(1,1);
		y4(i,1) = errors(i,1);
	else
		y1(i,1) = 0.6 * y1(i-1,1) + randn(1,1);
		y4(i,1) = errors(i,1) + 0.5 * errors(i-1,1);
	end;
	
	if i <=2 
		y2(i,1) = randn(1,1);
	else
		y2(i,1) = 1.3 * y2(i-1,1) - 0.5 * y2(i-2,1) + randn(1,1);
	end;
	
	if i <= 3
		y3(i,1) = randn(1,1);
		y5(i,1) = errors(i,1);

	else
		y3(i,1) = 1.3 * y3(i-1,1) - 0.75 * y3(i-2,1) + 0.3 * y3(i-3 ,1) + randn(1,1);
		y5(i,1) = errors(i,1) - 0.6 * errors(i-1,1) + 0.7 * y5(i-1 ,1) - 0.4 * y5(i-2 , 1);
	end;
	
	if i <= 4
		y6(i,1) = errors(i,1);
	else
		y6(i,1) = errors(i,1) + 0.7 * errors(i-1,1) - 0.5 * errors(i-2 ,1) + 0.3 * errors(i-3 , 1);
	end;
	
	
end;


save ar1.txt y1 -ASCII;
save ar2.txt y2 -ASCII;
save ar3.txt y3 -ASCII;
save ma1.txt y4 -ASCII;
save ma3.txt y6 -ASCII;
save arma.txt y5 -ASCII;

