%regression m-file using a design similar 
%to that appearing in the paper. 
y = 0;
clear all;
load float_all.txt;
%load c:\jia\brookings\finaldata\ptime.txt;
data=float_all;
date=data(:,1);
x=data(:,2);
ts=data(:,3);
fri=data(:,4);
dts1=data(:,5);
z1=data(:,6);
y=ts;
std=std(x);
mean=mean(x);
t1=prctile(x,25);
t2=prctile(x,75);
iqr=t2-t1;
a=min(std,(iqr/1.34));
h=0.9*a*(length(x))^(-0.2);
hmean = 0.8;
q = 0.5;
temp = ( 1/normcdf(q, 0, 1) )^2;
temp = normpdf(temp, 0, 1);
temp = (q * (1-q) / temp)^0.2;
hmed = hmean * temp;
points=linspace(4, 10, 601)';
bloclin = qreg(points,x,y,hmed);
plot(x,ts,'.');	
hold on;
%plot(x1,ts5,'*');
%hold on;
plot(bloclin(:,1),bloclin(:,2),'-');
%plot(points,frits,'-.');
xlabel('Time of Passing Toll Sign');	
ylabel('Time Saving');
legend('Data Points of Weekdays','Fitted Curve of Weekday');
hold off;
save time_50.txt bloclin -ASCII;

if y == 1
	%frits=bloclin;
	%for i=1:length(points);
   %if points(i)>=5.75
   	%frits(i)=frits(i)-1.413889;
   %end;
	%end;
	%clear i;

	%ts5=ts.*fri;
	%x1=zeros(length(x),1);
	%for i=1:length(x);
   	%if ts5(i)~=0
      	%x1(i)=x(i);
   	%end;
	%end;
	%ts5=nonzeros(ts5);
	%x1=nonzeros(x1);

	plot(x,ts,'.');	
	hold on;
	%plot(x1,ts5,'*');
	%hold on;
	plot(points,bloclin,'-');
	hold on;
	%plot(points,frits,'-.');
	xlabel('Time of Passing Toll Sign');	
	ylabel('Time Saving');
	legend('Data Points of Weekdays',...
   	'Fitted Curve of Weekday');
	%print -deps -tiff fig1;
	hold off;
	resid=y-bloclin;
	%plot(x,resid,'.');
	y=resid.*resid;
	points = linspace(4,10,601)';
	sigma=loclinear2(points,x,y,h);
	%y=ts+1.413889*z1;
	y=ts;
	bloclin=loclinear2(points,x,y,h);
	mean=[points,bloclin];
	variance=[points,sigma];
	%save mean.txt mean -ASCII;
	%save variance.txt variance -ASCII;
	%plot(points,bloclin,'-.');
	%hold on;
	%plot(points,sigma,'-');
	%xlabel('Time of passing Toll Sign');
	%ylabel('Value of Residual');
	%legend('Mean Time Saving(Weekday)','Conditional Variance of Time Saving',2);
	%print -deps -tiff fig3;
	%hold off;
	list_mean=zeros(length(ptime),1);
	list_std=zeros(length(ptime),1);
	step=0.0100;
	for i=1:601;
		p=4.000+(i-1)*step;   
   	for j=1:length(ptime);
      	if abs(ptime(j)-p)<0.0001
         	list_mean(j)=bloclin(i);
         	list_std(j)=sigma(i);
      	end;
   	end;
	end;
	save c:\jia\temp\etsave.txt list_mean -ASCII;
	save c:\jia\temp\estd.txt list_std -ASCII;
end;

