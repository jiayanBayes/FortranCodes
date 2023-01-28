function [rmat] = hrandn(m,n)
% the return value, rmat, is a NG x(NNC*REPTITION) matrix
%global REPTITION NNC

global NNC REPTITION;

rmat=zeros(m,1);
prime=[2;3;5;7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;...
	79;83;89;97;101;103;107;109;113];
i=1;
while i<=NNC,
	tt1=halton((m*REPTITION+10),prime(i));
	tt1=tt1(11:(m*REPTITION+10));
	x=zeros(m,1);
	j=1;
	while j<=REPTITION,
		gstart=(j-1)*m+1;
		gend=gstart+m-1;
		x=[x,tt1(gstart:gend)];
		j=j+1;
		clear gstart gend;
	end;
	x=x(:,2:(REPTITION+1));
	rmat=[rmat,x];
	i=i+1;
end;

rmat=rmat(:,2:(n+1));
        
	
			
