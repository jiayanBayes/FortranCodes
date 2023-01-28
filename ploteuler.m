
w = 1/0.995;
nsize = 500;

cvec = rand(nsize,1);
cvec = sort(cvec);

counter = 0;
for i = 1:nsize;
	counter = counter + 1;
	c = cvec(i,1);
	for j = 1:nsize;
		p = c + 0.1;
		
		if p > 1
			break;
		end;
		
		y = w - p + 1 - 4 * sqrt(w - p + 4 * (p - c) );
	end;
	
	temp = [c, p, y];
	
	if counter == 1
		res = temp;
	else
		res = [res;temp];
	end;
end;



