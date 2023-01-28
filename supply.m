function [t] = supply(c, v)

ft = 60/65;
l = 10;
a = 0.15;
b = 4;
ratio = v / c;
t = ft * l * (1 + a * ratio^b);

