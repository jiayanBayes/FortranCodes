clear all;
qupset = [10;20;30;40;50;60;70;80;90;100;110;120;130;140;150;160];
load ystar_qup_r08.txt;
load ystar_qup_rho10_r08.txt;
load ystar_qup_rho20_r08.txt;

%load ystar_qup_r10.txt;
%load ystar_qup_r12.txt;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 8.0 5]);

plot(qupset, ystar_qup_rho10_r08, 'kx-', 'LineWidth', 0.2);
hold on;
plot(qupset, ystar_qup_r08, 'ko-', 'LineWidth', 0.2);
hold on;
plot(qupset, ystar_qup_rho20_r08, 'k*-', 'LineWidth', 0.2);
hold off;
legend('{\rho} = 0.10','{\rho} = 0.15','{\rho} = 0.20');
ylabel('{\ity}_s^*', 'FontSize', 12);
xlabel('{\itnq}^+', 'FontSize', 12);

