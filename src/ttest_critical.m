function tcritical = ttest_critical(alpha,mu0,ndof)
% find the critical value for an F-test of equal variance
% 20200624 Kurt Feigl


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
% pN2 = cdf('norm',sqrt(variance_ratio),0,1);



% T-distribution
x=0.001:0.05:5;
Tc=cdf('T',x,ndof);
Tp=pdf('T',x,ndof);

% get bounds of T-statistic
%Flb = icdf('F',pNlb,ndof1,ndof2)
Tub = icdf('T',alpha,ndof)
tcritical = Tub;


pTmu0 = cdf('T',mu0,ndof);
Tmu0 = icdf('T',pTmu0, ndof)

confidence = 100*(1.0 - alpha);

% right-tail test
if pTmu0 > alpha
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
% plot([sqrt(variance_ratio), sqrt(variance_ratio)],[0, pN2],'c-');
xlabel('Z');
ylabel('P');
title('Standard normal (0, 1)');
legend('CDF','PDF','critical','test');

subplot(2,1,2);
hold on;
plot(x,Tc,'r-');
plot(x,Tp,'b-');
% %plot([Flb, Flb],[0, pNlb],'g-');
plot([Tub, Tub],[0, pNub],'g-');
plot([Tmu0, Tmu0],[0, pTmu0],'c-');
xlabel('x');
ylabel('P');
title(sprintf('T(%d) %s',ndof,test_string));
legend('CDF','PDF','critical','test');

print(gcf,'-dpdf',sprintf('%s.pdf',mfilename),'-r600','-fillpage','-painters');

return
end
