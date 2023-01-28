
clear all;
load Solutions_Bargaining_Revised_66.txt;
load Solutions_Bargaining_Revised_48.txt;
load Solutions_Bargaining_Revised_210.txt;

data1 = Solutions_Bargaining_Revised_66;
data2 = Solutions_Bargaining_Revised_48;
data3 = Solutions_Bargaining_Revised_210;

plot(data1(:,4), data1(:,3), 'r-');
hold on;
plot(data2(:,4), data2(:,3), 'b-');
hold on;
plot(data3(:,4), data3(:,3), 'g-');
hold off;
xlabel('Consumer surplus change relative to the no-toll base ($/person)', 'FontSize', 12);
ylabel('Operating profits of the monopoly ($/person)', 'FontSize', 12);
legend('{\itK}_{r1}=6000, {\itK}_{r2}=6000', '{\itK}_{r1}=4000, {\itK}_{r2}=8000', '{\itK}_{r1}=2000, {\itK}_{r2}=10000');

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
