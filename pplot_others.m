
%%% This code is for plotting the cost function under different parameter values %%%
clear all;
rset = [0.10;0.10;0.10];
rhoset = [0.15;0.15;0.15];
c = 50;
lamdset = [10;10;10];
sigmaset = [100;100;100];
qupset = [30;60;100];
qlow = 0;
kset = [1000;1000;1000];
ystar = zeros(length(kset), 1);

data = zeros(1000, length(kset) + 1); 
for i = 1:length(kset);
	r = rset(i,1);
	k = kset(i,1);
	sigma = sigmaset(i,1);
	rho = rhoset(i,1);
	lamd = lamdset(i,1);
	qup = qupset(i,1);
	[out, yoptim, dval] = mmi(r, rho, c, lamd, sigma, k, qup, qlow, 1);
	ystar(i,1) = yoptim;
	if i == 1
		data = out(:,1:3);
	else
		data = [data, out(:,2:3)]; 
	end;
end;

pdata = zeros(length(data), 4);
pdata(:,1) = data(:,1);

group = 1;
for i = 1:length(kset);
	for j = 1:length(data);
		if pdata(j,1) <= ystar(i,1)
			pdata(j, i+1) = data(j, group + 1);
		else
			pdata(j, i+1) = data(j, group + 2);	
		end;
	end;
	group = group + 2;
end; 
save pdata_qup.txt pdata -ASCII;
save ystar_qup.txt ystar -ASCII;
stop
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 4.0 5]);

plot(pdata(:,1), pdata(:,2), 'k-', 'LineWidth', 0.2);
hold on;
plot(pdata(:,1), pdata(:,3), 'k--', 'LineWidth', 0.2);
hold on;
plot(pdata(:,1), pdata(:,4), 'k-.', 'LineWidth', 0.2);
hold off;
legend('Inventory Cost = 0.05','Inventory Cost = 0.15','Inventory Cost = 0.30');
%axis([0,100,50,500]);
xlabel('Initial Reserve');
ylabel('Expected Cost');
print -deps -tiff fig2.eps
