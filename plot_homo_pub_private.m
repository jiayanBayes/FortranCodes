
clear all;

load homo_private_04.asc;
load homo_private_4525.asc;
load homo_private_53514.asc;

data1 = [homo_private_04; homo_private_4525; homo_private_53514];

load homo_pub.asc;
load homo_pub_0815.asc;

data2 = [homo_pub; homo_pub_0815];
