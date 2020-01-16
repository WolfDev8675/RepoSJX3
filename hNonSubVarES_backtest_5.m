function [uh, meane] = hNonSubVarES_5(r0, iu, fu)

%r1=flipud(Sample);  %Flip array up to down

temp4 = [];

clear r y y1 h1 h2 a s parm F1 F2 F3 t1 z1 z2 temp t45 Var5 V5 ES5 pvalue1 t11 t12 pd k k1;
r = datasample(r0,4000); %Sampling With Replacement from r1
y = sort(r);
k=myfunc(y,iu,fu);
%k = u;
k1=prctile(y,[99 95 90]);

[f12,x12]=ecdf(y);
f12=f12(2:end);
x12=x12(2:end);
y(y<=k)=[];
parm = gpfit(y-k);
[h1, h2] = size (y);
s(1)=0;
a(1)=0;
for i = 2:h1
    clear F1 F2 F3
    F1=1-interp1(x12,f12,y(i),'nearest');
    F2=1-interp1(x12,f12,k,'nearest');
%    F1= 1-normcdf(y(i),t11,t12);
%    F2 = 1- normcdf(k,t11,t12);
    F3 = F1/F2;
    g = 1-gpcdf(y(i)-k,parm(1), parm(2), 0);
    s(i)=(g-F3)^2;
    a(i) = a(i-1) + s(i);
end
y1=a';
%checking turning Point
t1=max(h1-15, 2);
temp=[];
temp1=[];
for j=t1:-1:t1-min(400, t1-1)
    clear z1 z2 b;
    z1 = y(j:h1);
    z2 = y1(j:h1); 
    z3 = [ones(h1-j+1,1) z1];
    b = regress(z2,z3);
    temp=[temp;[j y(j) y1(j) b(2)]];      
end
t3=temp(:,4)/std(temp(:,4));
t4=t3.^2;
pvalue=(1-chi2cdf(t4,1));
t45=[temp t4 pvalue];
cutoff=min(find(t45(:,6)<0.05));

Var5=y(t1-cutoff+1);
Exp5=mean(y(find(y>Var5)));

k2(1)= mean(y(find(y>k1(1))));

k2(2)= mean(y(find(y>k1(2))));

k2(3)= mean(y(find(y>k1(3))));

pvalue1 = 1-interp1(x12,f12,Var5,'nearest');

% k1 is the VaR at 99 95 90 confidence level & k2 is the Expected Shortfall in
% that level, k is the threshold of the non-subjective VaR, pvalue1 is
% the probability level of the non-subjective VaR, Var5 is the non-subjective VaR
% and Exp5 is the Expected Shortfall

temp4 = [temp4; [k1 k2 k pvalue1 Var5 Exp5]];

meane = temp4;
%stde = std(temp4);
%[~,~, mES5,stdES5] = hfindVaRES(r0,1-meane(5), 0, 1);

uh = meane(4);




