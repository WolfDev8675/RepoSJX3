function z2 = myfunc(r1,x1,x2)
warning('off', 'stats:gpfit:EvalLimit');
[bootstat,bootsam] = bootstrp(20,[],r1); %Sampling With Replacement from r1
Z=sort(r1(bootsam));% bootsam is the indices and Z1(bootsam) is 100 sample drawn  
temp = [];
for k1 = x1:.1:x2
a5=0;
for i = 1:20
      y5=Z(:,i);
    [f12,x12]=ecdf(y5);                                                                                                              
     f12=f12(2:end);
     x12=x12(2:end);
     k2 = min(k1,(y5(end)-0.1));
     k1 = k2;
    y5(y5<=k1)=[];
    para = gpfit(y5-k1);
    [h7, h8] = size (y5>k1);
    h4=min(find(y5>k1));
    [h5, h6] = size (y5);
    s5=0;
    clear F1 F2 F3 g;
    for j5=h4:h5
        
          F1=1-interp1q(x12,f12,y5(j5));
         F2=1-interp1q(x12,f12,k1);
        F3 = F1/F2;
        g = 1-gpcdf(y5(j5)-k1,para(1), para(2), 0);
        s5=s5+(g-F3)^2;
     end
 a5=a5+(s5/h7);
end
temp = [temp; [k1 a5]];
end
[M,I]= min(temp(:,2));
z2=temp(I,1);
end




