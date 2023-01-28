function [obf]=lae(b);
load tempp3.txt;
load yuse.txt;
load weights.txt;
q=0.5;
xb=tempp3*b;
resid=yuse-xb;
w=zeros(length(resid),1);
for j=1:length(resid);
   if resid(j)>0
      w(j)=q;
   else
      w(j)=(1-q);
   end;
end;
tp1=abs(resid);
tp2=w.*tp1;
tp3=weights.*tp2;
obf=sum(tp3);
