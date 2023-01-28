
clear all;

load plotobj_data.asc;
data = plotobj_data;
z = zeros(20,20);
y = zeros(20,0);
group = 0;
for i = 1:20;
	z(i,:) = data(group+1:group+20, 3);
	y(i,1) = data(group+1,1);
	group = group + 20;
end;
x = data(1:20,2);
mesh(x, y, z)
