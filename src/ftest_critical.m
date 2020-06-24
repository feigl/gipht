function fcritical = ftest_critical(alpha,variance_ratio,ndof1,ndof2)
% find the critical value for an F-test of equal variance
% 20200624 Kurt Feigl

if variance_ratio < 1
    error('variance ratio must be greater than unity for a right-tailed test');
end

% normal distribution
z=-4:0.1:4;
Nc = cdf('norm',z,0,1);
Np = pdf('norm',z,0,1);
% get probability of +/- one standard deviation
% pNlb = cdf('norm',-1,0,1); 
% pNub = cdf('norm',+1,0,1);
% get probability of alpha
%pNlb = alpha/2.;
pNub = 1. - alpha/2.;
%pNub = 1. - alpha;
zNub = icdf('norm',pNub,0,1);
% conf = pNub - pNlb
% get probability of variance ratio
pN2 = cdf('norm',sqrt(variance_ratio),0,1);



% Fdistribution
x=0.001:0.05:5;
Fc=cdf('F',x,ndof1,ndof2);
Fp=pdf('F',x,ndof1,ndof2);

% get bounds of F-statistic
%Flb = icdf('F',pNlb,ndof1,ndof2)
Fub = icdf('F',1-alpha,ndof1,ndof2)


pFra = cdf('F',variance_ratio,ndof1,ndof2);
Fra = icdf('F',pFra, ndof1,ndof2)

confidence = 100*(1.0 - alpha);

% right-tail test
if pFra > alpha
    H = 1;
    test_string = sprintf('Null hypothesis rejected with %.0f %% confidence',confidence)
else
    H = 0;
    test_string = sprintf('Null hypothesis fails to be rejected with %.0f %% confidence',confidence);
end
   


figure; 
subplot(2,1,1);
hold on;
plot(z,Nc,'r-');
plot(z,Np,'b-');
%plot([-1 -1],[0, pNlb],'g-');
%plot([+1 +1],[0, pNub],'g-');
plot([zNub zNub],[0, pNub],'g-');
plot([sqrt(variance_ratio), sqrt(variance_ratio)],[0, pN2],'c-');
xlabel('Z');
ylabel('P');
title('Standard normal (0, 1)');
legend('CDF','PDF','critical','test');

subplot(2,1,2);
hold on;
plot(x,Fc,'r-');
plot(x,Fp,'b-');
% %plot([Flb, Flb],[0, pNlb],'g-');
plot([Fub, Fub],[0, pNub],'g-');
plot([Fra, Fra],[0, pFra],'c-');
xlabel('x');
ylabel('P');
title(sprintf('F(%d, %d) %s',ndof1, ndof2,test_string));
legend('CDF','PDF','critical','test');

fcritical = Fub;
print(gcf,'-dpdf',sprintf('%s.pdf',mfilename),'-r600','-fillpage','-painters');

return
end
