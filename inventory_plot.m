
w = 1/0.995;
p = 0.8;  
diff = 0.4;

yini = -0.5 * (w - p) + 0.5 * sqrt( (p-w)^2 + 4 * diff );
kstar = yini;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 8.0 5]);

%subplot(1,2,1);
xvec = linspace(0.001,1,1000)';
plot(xvec, xvec, 'k-');
hold on;

y = zeros(length(xvec),1);
for j = 1:length(xvec);
	x = xvec(j,1);
	y(j,1) = (w * diff / x) + p - w;
end; 

plot(xvec, y, 'k-');
hold on;

equi = yini * ones(length(xvec),1);
h1 = linspace(0, yini, length(xvec))';
plot(equi, h1, 'k:');
hold on;
plot(h1, equi, 'k:');
hold off;
axis([0,1,0,1]);

xlabel('{\itF}_d({\itK}_{\itt})');
ylabel('{\itF}_d({\itK}_{\itt+1})');
text(0.5, 0.5, '{\delta}', 'FontSize', 16);
text(0.2, 0.8, '{\itK}^*')
text(0.1, 0.5, '{', 'FontSize', 25);
text(0.8, 0.2, '{\itG}({\itF}_d({\itK}_t))');
text(0.2, 0.8, '{\it4}{\it5}^o');
text(0.1, 0.5, '(a)');

%subplot(1,2,2);
time = linspace(0, 20, 100)';
for iloop = 1:3;
	if iloop == 1
		w = 1/0.995;
	elseif iloop == 2
		w = 1/0.95;
	else
		w = 1/0.90;
	end;
	
	yini = -0.5 * (w - p) + 0.5 * sqrt( (p-w)^2 + 4 * diff );
	kstar = yini;

	nperiods = 2000;
	res = [0, yini];
	for i = 1: nperiods;
		yupdate = (w / yini) * diff + p - w;
		if yupdate > 1
			yupdate = 1;
		end;
	
		temp = [i, yupdate];
		res = [res; temp];
		yini = yupdate;
	
		if yupdate == 1
			break;
		end;
	end;

	if iloop == 1
		subplot(3,1,1);
	elseif iloop == 2
		subplot(3,1,2);
	else
		subplot(3,1,3);
	end;
	
	plot(res(:,1), res(:,2), 'k-');
	hold on;
	temp = kstar * ones(length(time),1);
	plot(time, temp, 'k:')
	hold off;
	axis([0, 20, 0, 1]);
	
	if iloop == 3
		xlabel('Time Period: {\itt}');
	end;
	ylabel('{\itF}_d({\itK}_t)');
	%text(0.2, 0.8, '{\itF}_d({\itK}^*)')
		
	%if iloop == 1
	%	text(0.9,0.1,'(a)');
	%elseif iloop == 2
	%	text(0.9,0.1,'(b)');
	%else
	%	text(0.9,0.1,'(c)');
	%end;
end;

%print -depsc -tiff -r600 fig1.eps
	