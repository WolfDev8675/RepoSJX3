r10 = xlsread('data.xlsx', 'loss', 'F2:F80000'); % Russia MICEX (F2:F80000 because this command will read only nonzero elements in the column)
loss6 = flipud(r10);  %Flip array up to down

T = size(loss6, 1)                         %sample size
T1=round(T/2);
T2=T-T1;

j8=0


for j6=[1,9,19,29,36,78]

    j6
loss5=loss6(j6:T1+j6-1);

% Model is fitted based on rolling windows of T1 observations

%  Finding ARMA order for best fit using AIC and BIC criteria
ai = zeros(4,4); % Initialize
bi = zeros(4,4);
% try-catch-continue clause is for going to next iteration when there is a
% FATAL ERROR
for p = 1:4
    for q = 1:4
        moda = arima('Constant', 0,'ARLags',  [1:p], 'MALags', [1:q], 'Distribution', 't');
        %options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',500, 'TolFun',.05, 'Display','Off');
        try
            [fita,~,logLa] = estimate(moda,loss5, 'Display','Off');
        catch
           disp(['error at j6=' j6]);
           continue
        end
        
        [aic, bic]= aicbic (logLa,p+q,T);
        ai(p,q) = aic;
        bi(p,q) = bic;
     end
end

[ia,ja] = find(ai == min(ai(:)));
[ia1,ja1] = find(bi == min(bi(:)));

%fprintf('AIC Values): %d\t %d\n', ia,ja)
%fprintf('BIC Values): %d\t %d\n', ia1,ja1)


% Computing ARMA REsidual

Mdl = arima ('Constant', 0,'ARLags',  [1:ia], 'MALags', [1:ja], 'Distribution', 't');
%options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',500, 'TolFun',.05, 'Display','Off');
try
    EstMdl = estimate(Mdl,loss5, 'Display','Off'); 
catch
     disp(['error at j6=' j6]);
     continue
end
        
        

% 'Infer' computes residual 

[E]=infer(EstMdl,loss5);

mut=loss5-E;   % mean of the ARMA model    



%Finding GARCH order for best fit using AIC criteria

for p = 1:4
    for q = 1:4
        if j6==1
            modg   = egarch('Constant', NaN,'GARCHLags', [1:p], 'ARCHLags',[1:q],'Distribution', 't');
        %modg = egarch(p,q);
        options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',100,'StepTolerance',0.03,'OptimalityTolerance',0.02, 'ConstraintTolerance',.02, 'Display','Off');
        try
            [fitg,~,logLg,ifn] = estimate(modg, E, 'Display','Off');
        catch EXC
           disp(['error at j6=' j6]);
           continue
        end
        else
        modg   = egarch('Constant', NaN,'GARCHLags', [1:p], 'ARCHLags',[1:q],'Distribution', 't');
        %modg = egarch(p,q);
        options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',100,'StepTolerance',0.03,'OptimalityTolerance',0.02, 'ConstraintTolerance',.02,'Display','Off');
        prevMod=[fitg.Constant,cell2mat(fitg.GARCH),cell2mat(fitg.ARCH)]';
        try
            [fitg,~,logLg] = estimate(modg, E,'E0',ifn.X ,'Display','Off');
        catch EXC
           disp(['error at j6=' j6]);
           continue
        end
        end
         [aic5, bic5]= aicbic (logLg,p+q,T);
        aig(p,q) = aic5;
        big(p,q) = bic5;
        
     end
end

[iag,jag] = find(aig == min(aig(:)));
[iag1,jag1] = find(big == min(big(:)));

%fprintf('AIC Values for Volatility model): %d\t %d\n', iag,jag)
%fprintf('BIC Values for Volatility model): %d\t %d\n', iag1,jag1)

%
if j6==1
    md2   = egarch('Constant', NaN,'GARCHLags', [1:iag], 'ARCHLags',[1:jag],'Distribution', 't');
    options1 = optimoptions(@fmincon,'Diagnostics','off','Algorithm','sqp','MaxIter',150,'StepTolerance',0.03,'OptimalityTolerance',0.02, 'ConstraintTolerance',.02, 'Display','Off');
    try
    [EstMd2,~,~,ifm] = estimate(md2, E, 'options', options1,'Display','Off');
    catch EXC
    disp(['error at j6=' j6]);
    continue
end
else
    md2   = egarch('Constant', NaN,'GARCHLags', [1:iag], 'ARCHLags',[1:jag],'Distribution', 't');
    options1 = optimoptions(@fmincon,'Diagnostics','off','Algorithm','sqp','MaxIter',100,'StepTolerance',0.03,'OptimalityTolerance',0.02, 'ConstraintTolerance',.02, 'Display','Off');
    prev=[EstMd2.Constant,cell2mat(EstMd2.GARCH),cell2mat(EstMd2.ARCH)]';
    try
    [EstMd2,~,~,ifm] = estimate(md2, E,'E0',ifm.X,'V0',V, 'options', options1,'Display','Off');
    catch EXC
    disp(['error at j6=' j6]);
    continue
end
end

%

        
        
%options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',200, 'TolFun',.05, 'Display','Off');
%EstMd2 = estimate (md2, E, 'Display','Off');

% Estimation of Conditional Variance

[V]=infer(EstMd2,E);

scg=sqrt(V); %standard DEviation

Error5 = E./scg;

% Sampling With Replacement - Difference between Ordinary ARMA-GARCH model
% and FHS model

E6 = datasample(Error5,10000); %Sampling With Replacement from ERROR5


%Initial estimate of the Threshold (iu); Minimum 500 data point is required
%for fitting GPD



y90 = sort(E6);

g89= y90(end-500);

iu = myfunc(y90,0,g89);




% Choosing a small neighbourhood of iu (-0.2 to +0.2) for computing the
% Non-subjective Value-at-Risk
iuz = iu - 0.2;
fuz = iu + 0.2;
[u, meana] = hNonSubVarES_backtest_5(y90, iuz, fuz);

% temp4: k1 - ConvenTional VaR(99,95,90); k2 - Expected Shortfall(99,95,90);
% k - threshold; pvalue1 - probability; Var5 -Non-Subjective Value-at-Risk;
% Exp5 - non-Subjective Expected Shortfall

meanf = mut(end) + (scg(end))*meana;

mVarhNS = mut(end) + (scg(end))*meana(9);
%stdVarhNS = sigt*stda(6);
mEShNS = mut(end) + (scg(end))*meana(10);
%stdEShNS = sigt*stda(7);

% V1 is the non-subjective VaR forecast
% V2 is next days loss
V1(j6)=mVarhNS;
V2(j6)=loss6(T1+j6);

 if (V2(j6)>=V1(j6))
       j8=j6+j8; %Total Number of violation
 end
   
end

j8
T2
   
j8/T2  %proportion of violation

%prob1= 1-poisscdf(j6,sum1)




