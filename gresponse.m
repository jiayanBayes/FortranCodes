
clear all;

load profitresponse_equalcapacity_part1.asc;
load profitresponse_equalcapacity_part2.asc;
data1 = profitresponse_equalcapacity_part1;
data2 = profitresponse_equalcapacity_part2;

%subplot(3,1,1);
plot(data1(:,2), data1(:,1), 'r-');
hold on;
plot(data1(:,1), data1(:,2), 'b-');
hold on;
plot(data2(:,2), data2(:,1), 'r-');
hold on;
plot(data2(:,1), data2(:,2), 'b-');
hold off;
ylabel('Toll on Route r1($)', 'FontSize', 12);
xlabel('Toll on Route r2($)', 'FontSize', 12);
text(10,10, 'Equilibrium', 'FontSize', 10);
text(8,10, 'Equilibrium', 'FontSize', 10);
legend('Response Curve of Operator r2', 'Response Curve of Operator r1', 'Location', 'SouthWest');
title('Solutions of Betrand competition');

%load profitfunction_leader_duopoly_part1.asc;
%load profitfunction_leader_duopoly_part2.asc;
%data1 = profitfunction_leader_duopoly_part1;;
%data2 = profitfunction_leader_duopoly_part2;
%subplot(2,1,2);
%plot(data1(1:108,1), data1(1:108,2), 'r');
%hold on;
%plot(data2(:,1), data2(:,2), 'r');
%hold off;
%xlabel('Toll on Route r2 ($)', 'FontSize', 12);
%ylabel('Profit of Operator r2 as the Price Leader ($)', 'FontSize', 11);
%title('b: Profit function of the price leader in Stackelberg competition');


load response_travelersolution_part1.asc;
load response_travelersolution_part2.asc;
data1 = response_travelersolution_part1;
data2 = response_travelersolution_part2;

%subplot(2,1,1);
plot(data1(:,2), data1(:,1), 'r-');
hold on;
plot(data1(:,1), data1(:,2), 'b-');
hold on;
plot(data2(:,2), data2(:,1), 'r-');
hold on;
plot(data2(:,1), data2(:,2), 'b-');
hold off;
ylabel('Toll on Route r1($)', 'FontSize', 12);
xlabel('Toll on Route r2($)', 'FontSize', 12);
text(3,3, 'Equilibrium', 'FontSize', 10);
text(5,3, 'Equilibrium', 'FontSize', 10);
legend('Response Curve of Operator r1', 'Response Curve of Operator r2', 'Location', 'NorthEast');
axis([0,10,0,10]);
title('Solutions to Bertrand competition with bargaining: traveler solution');
%stop

%load welfare_leader_revision_privatebargaining_part1.asc;
%load welfare_leader_revision_privatebargaining_part2.asc;
%data1 = welfare_leader_revision_privatebargaining_part1;
%data2 = welfare_leader_revision_privatebargaining_part2;
%subplot(2,1,2);
%plot(data1(:,1), data1(:,2)/24300 - 12.12, 'r');
%hold on;
%plot(data2(:,1), data2(:,2)/24300 - 12.12, 'r');
%hold off;
%xlabel('Toll on Route r2 ($)', 'FontSize', 12);
%ylabel('Consumer Surplus Change ($/person)', 'FontSize', 12);
%title('b: Consumer surplus w.r.t. the toll of the price leader under Bargaining Stackelberg competition: traveler solution');
%axis([0,10, -0.5, 1]);

load Response_Bargaining_Profit_Part1.txt;
load Response_Bargaining_Profit_Part2.txt;
data1 = Response_Bargaining_Profit_Part1;
data2 = Response_Bargaining_Profit_Part2;

subplot(1,1,1);
plot(data1(:,1), data1(:,2), 'b');
hold on;
plot(data1(:,2), data1(:,1), 'r');
hold on;
plot(data2(:,1), data2(:,2), 'b');
hold on;
plot(data2(:,2), data2(:,1), 'r');
hold off;
ylabel('Toll on Route r1($)', 'FontSize', 12);
xlabel('Toll on Route r2($)', 'FontSize', 12);
%text(3,3, 'Equilibrium A', 'FontSize', 10);
%text(5,3, 'Equilibrium B', 'FontSize', 10);
legend('Response Curve of Operator r2', 'Response Curve of Operator r1', 'Location', 'NorthEast');
title('Multiple equilibria of Bertrand competition with bargaining: operator solution');
%axis([11.688, 11.6892, 2.231 2.2318]);

%load profitleader_bargaining_part1.asc;
%load profitleader_bargaining_part2.asc;
%data1 = profitleader_bargaining_part1;
%data2 = profitleader_bargaining_part2;
%subplot(3,1,2);
%plot(data1(:,1), data1(:,3)/24300, 'r');
%hold on;
%plot(data2(:,1), data2(:,3)/24300, 'r');
%hold off;
%xlabel('Toll on Route r2 ($)', 'FontSize', 12);
%ylabel('Profit of Operator r2 ($/person)', 'FontSize', 10);
%title('b: Profit of price leader in bargaining Stackelberg competition: operator solution');

%subplot(3,1,3);
%plot(data1(:,1), data1(:,2)/24300 - 12.12, 'r');
%hold on;
%plot(data2(:,1), data2(:,2)/24300 - 12.12, 'r');
%hold off;
%xlabel('Toll on Route r2 ($)', 'FontSize', 12);
%ylabel('Consumer Surplus Change ($/person)', 'FontSize', 10);
%title('c: Consumer surplus w.r.t. the toll of the price leader under Bargaining Stackelberg competition: operator solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load response_operatorsolution_48_part1.asc;
%load response_operatorsolution_48_part2.asc;
%data1 = response_operatorsolution_48_part1;
%data2 = response_operatorsolution_48_part2;

%load response_operatorsolution_84_part1.asc;
%load response_operatorsolution_84_part2.asc;
%data3 = response_operatorsolution_84_part1;
%data4 = response_operatorsolution_84_part2;

%subplot(1,1,1);
%plot(data1(:,2), data1(:,1), 'r');
%hold on;
%plot(data4(:,2), data4(:,1), 'b', 'LineWidth', 2);
%stop
%subplot(1,1,1);
%plot(data1(:,1), data1(:,2), 'r');
%hold on;
%plot(data3(:,2), data3(:,1), 'b', 'LineWidth', 2);
%hold on;
%plot(data2(:,1), data2(:,2), 'r');
%hold on;
%plot(data4(:,2), data4(:,1), 'b', 'LineWidth', 2);
%hold off;
%ylabel('Toll on Route r2($)', 'FontSize', 12);
%xlabel('Toll on Route r1($)', 'FontSize', 12);
%axis([2.55, 2.65, 9.2, 9.45]);

%stop
%text(3,3, 'Equilibrium (this equilibrium constitutes the unique SPE to the overall game)', 'FontSize', 10);
%text(5,3, 'Equilibrium', 'FontSize', 10);
%legend('Response Curve of Public Operator (Operator r2)', 'Response Curve of Private Operator (Operator r1)', 'Location', 'NorthWest');
%title('a: Solution for Bertrand competition');

%stop
%%%%%%%%%%%%%%%%%%%%%%%%%
load pubresponse_capacity102_part1.asc;
load pubresponse_capacity102_part2.asc;
data1 = pubresponse_capacity102_part1;
data2 = pubresponse_capacity102_part2;

load profitresponse_capacity210_part1.asc;
load profitresponse_capacity210_part2.asc;
data3 = profitresponse_capacity210_part1;
data4 = profitresponse_capacity210_part2;

subplot(1,1,1);
plot(data1(:,1), data1(:,2), 'r');
hold on;
plot(data3(:,2), data3(:,1), 'b', 'LineWidth', 2);
hold on;
plot(data2(:,1), data2(:,2), 'r');
hold on;
plot(data4(:,2), data4(:,1), 'b', 'LineWidth', 2);
hold off;
ylabel('Toll on Route r2($)', 'FontSize', 12);
xlabel('Toll on Route r1($)', 'FontSize', 12);
text(3,3, 'Equilibrium A (constitutes the unique SPE to the overall game)', 'FontSize', 10);
text(5,3, 'Equilibrium B', 'FontSize', 10);
legend('Response Curve of Public Operator (Operator r2)', 'Response Curve of Private Operator (Operator r1)', 'Location', 'NorthWest');
title('Solution for Bertrand public-private competition');

load response_travelersolution_48_part1.asc;
load response_travelersolution_48_part2.asc;
data1 = response_travelersolution_48_part1;
data2 = response_travelersolution_48_part2;

load response_travelersolution_84_part1.asc;
load response_travelersolution_84_part2.asc;
data3 = response_travelersolution_84_part1;
data4 = response_travelersolution_84_part2;

subplot(1,1,1);
plot(data1(:,2), data1(:,1), 'r');
hold on;
plot(data3(:,1), data3(:,2), 'b');
hold on;
plot(data2(:,2), data2(:,1), 'r');
hold on;
plot(data4(:,1), data4(:,2), 'b');
hold off;
ylabel('Toll on Route r1($)', 'FontSize', 12);
xlabel('Toll on Route r2($)', 'FontSize', 12);
text(3,3, 'Equilibrium', 'FontSize', 10);
text(5,3, 'Equilibrium', 'FontSize', 10);
legend('Response Curve of Operator r1', 'Response Curve of Operator r2', 'Location', 'NorthEast');
title('Solution for Bertrand bargaining competition with unequal capacity allocation ({\itK}_{\itr1}=4000,{\itK}_{\itr2}=8000): traveler solution');
stop

load response_operatorsolution_48_part1.asc;
load response_operatorsolution_48_part2.asc;
data1 = response_operatorsolution_48_part1;
data2 = response_operatorsolution_48_part2;

load response_operatorsolution_84_part1.asc;
load response_operatorsolution_84_part2.asc;
data3 = response_operatorsolution_84_part1;
data4 = response_operatorsolution_84_part2;

subplot(1,1,1);
plot(data1(:,2), data1(:,1), 'r');
hold on;
plot(data3(:,1), data3(:,2), 'b');
hold on;
plot(data2(:,2), data2(:,1), 'r');
hold on;
plot(data4(:,1), data4(:,2), 'b');
hold off;
ylabel('Toll on Route r1($)', 'FontSize', 12);
xlabel('Toll on Route r2($)', 'FontSize', 12);
%text(3,3, 'Equilibrium A (constitutes the unique SPE to the overall game)', 'FontSize', 10);
%text(5,3, 'Equilibrium B', 'FontSize', 10);
legend('Response Curve of Operator r1', 'Response Curve of Operator r2', 'Location', 'NorthEast');
title('Multiple equilibria for Bertrand bargaining competition with unequal capacity ({\itK}_{\itr1}=4000,{\itK}_{\itr2}=8000): operator solution');

%load pubfunction_leader_duopoly_part1.asc;
%load pubfunction_leader_duopoly_part2.asc;
%data1 = pubfunction_leader_duopoly_part1;
%data2 = pubfunction_leader_duopoly_part2;
%subplot(2,1,2);
%plot(data1(:,1), data1(:,2), 'r');
%hold on
%plot(data2(:,1), data2(:,2), 'r');
%hold off;
%xlabel('Toll on Route r2 ($)', 'FontSize', 12);
%ylabel('Welfare of Public Operator (Operator r2) as the Price Leader ($)', 'FontSize', 11);
%title('b: Welfare function of the public operator in Stackelberg competition');
