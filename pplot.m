
%%% This code is for plotting the cost function under different parameter values %%%
clear all;
r = 0.10;
rho = 0.15;
c = 50;
lamd = 10;
sigma = 100;
qup = 60;
qlow = 0;
Kset = [1000];
ystar = zeros(length(Kset), 1);

data = zeros(1000, length(Kset) + 1); 
for i = 1:length(Kset);
	k = Kset(i,1);
	[out, yoptim, dval] = mmi(r, rho, c, lamd, sigma, k, qup, qlow, 1);
	ystar(i,1) = yoptim;
	if i == 1
		data = out;
	else
		data = [data, out(:,2:3)]; 
	end;
end;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 8.0 5]);

star = yoptim * ones(length(data),1);
opti = ecost(yoptim) * ones(length(data),1);
temp1 = linspace(0, ecost(yoptim), length(data))';
temp2 = linspace(0, yoptim, length(data))';

%plot(data(:,1), data(:,2), '-');
%xlabel('Initial Reserve');
%ylabel('Expected Cost');
%axis([0,100,min(data(:,2)),1200]);
%print -deps -tiff -r300 fig1a.eps;

%plot(data(:,1), data(:,3), '--');
%xlabel('Initial Reserve');
%ylabel('Expected Cost');
%print -deps -tiff -r300 fig1b.eps;

plot(data(:,1), data(:,2), 'k-', 'LineWidth', 0.25);
hold on;
plot(data(:,1), data(:,3), 'k--','LineWidth', 0.25);
hold on;
plot(star, temp1, 'k:','LineWidth', 0.25);
hold on;
plot(temp2, opti, 'k:','LineWidth', 0.25);
hold off;
xlabel('{\ity}_s');
ylabel('{\itV( {\ity}_s )}');
axis([0,150,400,1200]);
text(85, 800, '{\ity}_s^* = 77.84, {\itV({\ity}_s^*)} = 718.49', 'EdgeColor', 'k');
text(25, 1000, '{\itV}^{\itL}({\ity}_s)');
text(125, 550, '{\itV}^{\itU}({\ity}_s)');
%print -depsc -tiff -r600 fig1.eps

