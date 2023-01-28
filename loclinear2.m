%This fctn computes local linear regression estimates using the BIWEIGHT kernel. 
%It applies to the regession model:
%
%                       y = m(x) + u,           where x is a scalar variable.
%
%The inputs to the m-file are:
%               the set of points where the function will be estimated (points)
%               the x variable  
%               the y variable
%               the bandwidth (h)
%NOTE THAT POINTS, x,and y SHOULD BE SPECIFIED AS COLUMN VECTORS.   
%
%The syntax is
%function out = loclinear2(points,x,y,h);
%
%The choice of bandwidth is not automatic, but is left to the 
%user to specify. See Fan and Gijbels (1996) Local Polynomial 
%Modelling and its Applications, Chapter 4 for alternate bandwidth 
%selection procedures.
%
%The output (out) is the CMF. An estimate of the CMF is obtained for 
%each value in the points vector.

function out = loclinear2(points,x,y,h);
n = length(points);
const = ones(length(x),1);
bhat = zeros(2,length(points));
predict = zeros(n,1);
for i = 1:n;
   a = -.5*sign( abs( (points(i,1)*const - x)/h ) -1 ) + .5;    %get the right data points, (K(x) ~=0)
   tempp0=nonzeros(x.*a);
   tempp1=nonzeros(y.*a);
   xandy =[tempp0,tempp1];                                                   
   ztheta = xandy(:,1);   
   yuse = xandy(:,2);
   q1 = (ztheta - points(i,1));
   nt2 = ( (ztheta- points(i,1))/h );
   weights = diag( (15/16)* ( 1-(nt2.^2)).^2 );
   tempp2=ones(length(ztheta),1);     
   tempp3 = [tempp2, q1];
   bhat(:,i) = inv(tempp3'*weights*tempp3)*tempp3'*weights*yuse;
   clear xandy ztheta yuse;
end;
out = bhat(1,:)';

