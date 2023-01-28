load ystar_qup.txt;
load pdata_qup.txt;
ystar_K = ystar_qup;
pdata_K = pdata_qup; 

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 8.0 5]);

star1 = ystar_K(1,1) * ones(length(pdata_K),1);
star2 = ystar_K(2,1) * ones(length(pdata_K),1);
star3 = ystar_K(3,1) * ones(length(pdata_K),1);
temp1 = linspace(0, 800.27, length(pdata_K))';
temp2 = linspace(0, 718.49, length(pdata_K))';
temp3 = linspace(0, 621.11, length(pdata_K))';

opti1 = 800.27 * ones(length(pdata_K),1);
h1 = linspace(0, ystar_K(1,1), length(pdata_K))';
opti2 = 718.49 * ones(length(pdata_K),1);
h2 = linspace(0, ystar_K(2,1), length(pdata_K))';
opti3 = 621.11 * ones(length(pdata_K),1);
h3 = linspace(0, ystar_K(3,1), length(pdata_K))';

plot(pdata_K(:,1), pdata_K(:,2), 'k-', 'LineWidth', 0.2);
hold on;
plot(pdata_K(:,1), pdata_K(:,3), 'k--', 'LineWidth', 0.2);
hold on;
plot(pdata_K(:,1), pdata_K(:,4), 'k-.', 'LineWidth', 0.2);
hold on;
plot(star1, temp1, 'k:', 'LineWidth', 0.2);
hold on;
plot(h1, opti1, 'k:', 'LineWidth', 0.2);
hold on;
plot(star2, temp2, 'k:', 'LineWidth', 0.2);
hold on;
plot(h2, opti2, 'k:', 'LineWidth', 0.2);
hold on;
plot(star3, temp3, 'k:', 'LineWidth', 0.2);
hold on;
plot(h3, opti3, 'k:', 'LineWidth', 0.2);
hold off;
legend('{\itq}^+ = 30','{\itq}^+ = 60','{\itq}^+ = 100');
xlabel('{\ity}_s');
ylabel('{\itV( {\ity}_s )}');
axis([0, 300, 300, 1000]);
text(25,  600, '{\ity}_s^* = 132.47, {\itV({\ity}_s^*)} = 800.27', 'EdgeColor', 'k');
text(100, 600, '{\ity}_s^* = 77.84, {\itV({\ity}_s^*)} = 718.49', 'EdgeColor', 'k');
text(50,  900, '{\ity}_s^* = 12.73, {\itV({\ity}_s^*)} = 621.11', 'EdgeColor', 'k');

%print -deps -tiff fig2.eps
