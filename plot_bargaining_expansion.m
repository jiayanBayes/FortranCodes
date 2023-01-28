
clear all;
load Solutions_Bargaining_Revised_86.txt;
load Solutions_Bargaining_Revised_410.txt;
load Solutions_Bargaining_Revised_212.txt;

data1 = Solutions_Bargaining_Revised_86;
data2 = Solutions_Bargaining_Revised_410;
data3 = Solutions_Bargaining_Revised_212;

plot(data1(:,4), data1(:,3), 'r-');
hold on;
plot(data2(:,4), data2(:,3), 'b-');
hold on;
plot(data3(:,4), data3(:,3), 'g-');
hold off;
xlabel('Consumer surplus change relative to the no-toll no-expansion base ($/person)', 'FontSize', 12);
ylabel('Operating profits of the monopoly ($/person)', 'FontSize', 12);
legend('{\itK}_{r1}=6000, {\itK}_{r2}=8000', '{\itK}_{r1}=4000, {\itK}_{r2}=10000', '{\itK}_{r1}=2000, {\itK}_{r2}=12000');
axis([0,3.11,0,3.5]);

x = data3(:,4);
y = data3(:,3);
ix0 = linspace(0, 3.1, 50)';
iy0 = zeros(50,1);

ix07 = linspace(0.7, 3.1, 50)';
iy07 = zeros(50,1);

ix14 = linspace(1.40, 3.1, 50)';
iy14 = zeros(50,1);

for i =1:50;
	iy0(i) = interp1(x, y, ix0(i));
	iy0(i) = iy0(i) - 3.18 - 0.12;

	iy07(i) = interp1(x, y, ix07(i));
	iy07(i) = iy07(i) - 2.71 - 0.12;

	iy14(i) = interp1(x, y, ix14(i));
	iy14(i) = iy14(i) - 0.89 - 0.12;
end;

plot(ix0, iy0, 'r');
hold on;
plot(ix07, iy07, 'b');
hold on;
plot(ix14, iy14, 'g');
hold off;
xlabel('Consumer surplus change relative to the no-toll no-expansion base ($/person)', 'FontSize', 12);
ylabel('Marginal profits of adding one more lane ($/person)', 'FontSize', 12);
legend('Initial outcome: operator solution', 'initial outcome: intermediate solution', 'initial outcome: traveler solution');
axis([0,3.11,0,2.5]);
stop

plot(data3(:,4), data3(:,1), 'r-', 'LineWidth', 2);
hold on;
plot(data3(:,4), data3(:,2), 'b-');
hold off;
legend('toll on route r1 ({\itK}_{r1}=2000)', 'toll on route r2 ({\itK}_{r2}=10000)');
xlabel('Consumer surplus change relative to the no-toll base ($/person)', 'FontSize', 12);
ylabel('Toll ($)', 'FontSize', 12);
axis([0,1.4,0,15]);
hold off;

x = data3(:,3) * 24300 * 250/0.045;
x = x / 10000000
plot(data3(:,4), x, 'r-')
xlabel('Consumer surplus change relative to the no-toll base ($/person)', 'FontSize', 12);
ylabel('Sell pricel ($ million per mile)', 'FontSize', 12);

stop

subplot(3,1,1)
plot(data1(:,4), data1(:,1), 'r-', 'LineWidth', 10);
hold on;
plot(data1(:,4), data1(:,2), 'b-');
hold off;
legend('toll on route r1', 'toll on route r2');
xlabel('Consumer surplus change relative to the no-toll base ($/person)');
ylabel('$')
axis([0,1.4,0,15]);
title('(a). Capacity allocation: {\itK}_{r1}=6000, {\itK}_{r2}=6000 ');

subplot(3,1,2)
plot(data2(:,4), data2(:,1), 'r-');
hold on;
plot(data2(:,4), data2(:,2), 'b-');
hold off;
legend('toll on route r1', 'toll on route r2');
xlabel('Consumer surplus change relative to the no-toll base ($/person)');
ylabel('$')
axis([0,1.4,0,15]);
title('(b). Capacity allocation: {\itK}_{r1}=4000, {\itK}_{r2}=8000 ');

subplot(3,1,3)
plot(data3(:,4), data3(:,1), 'r-');
hold on;
plot(data3(:,4), data3(:,2), 'b-');
hold off;
legend('toll on route r1', 'toll on route r2');
xlabel('Consumer surplus change relative to the no-toll base ($/person)');
ylabel('$')
axis([0,1.4,0,15]);
title('(c). Capacity allocation: {\itK}_{r1}=2000, {\itK}_{r2}=10000 ');
hold off;
