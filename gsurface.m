
clear all;
load plotobj_data_profit_euqalcapacity.asc;
data1 = plotobj_data_profit_euqalcapacity;

load plotobj_data_wel_euqalcapacity.asc;
data2 = plotobj_data_wel_euqalcapacity;

load plotobj_data_profit.asc;
data3 = plotobj_data_profit;

load plotobj_data_wel.asc;
data4 = plotobj_data_wel;

x1 = linspace(0, 30, 20)'; 
x2 = linspace(0, 30, 20)';

[i,j] = size(data1);
cons = 10 * ones(i, j);
mesh(x1, x2, data1); 
hold on;
mesh(x1, x2, cons);
hold off;
