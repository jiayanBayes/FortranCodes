function out = qreg(points,x,y,h);
n = length(points);
const = ones(length(x),1);
out = zeros(n,2);
for i = 1:n;
   a = -.5*sign( abs( (points(i,1)*const - x)/h ) -1 ) + .5;    %get the right data points, (K(x) ~=0)
   tempp0=nonzeros(x.*a);
   tempp1=nonzeros(y.*a);
   xandy =[tempp0,tempp1];                                                   
   ztheta = xandy(:,1);   
   yuse = xandy(:,2);
   q1 = (ztheta - points(i,1));
   nt2 = ( (ztheta- points(i,1))/h );
   weights = (15/16)* ( 1-(nt2.^2)).^2;
   tempp2=ones(length(ztheta),1);     
   tempp3 = [tempp2, q1];
   b0 = inv(tempp3'*tempp3)*tempp3'*yuse;
   save tempp3.txt tempp3 -ASCII;
   save yuse.txt yuse -ASCII;
   save weights.txt weights -ASCII;
   %save bo.txt b0 -ASCII;
   options=optimset('MaxFunEvals',300000000);                                             
   bhat=fminsearch('lae', b0, options);
   save bhat.txt bhat -ASCII;
   out(i,:)= [points(i,1), bhat(1)];
   clear xandy ztheta yuse b xb resid obf w tp1 tp2 bhat;
end;

