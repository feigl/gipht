
% normal distribution
z=-4:0.1:4;
Nc = cdf('norm',z,0,1);
Np = pdf('norm',z,0,1);
% get probability of +/- one standard deviation
pNlb = cdf('norm',-1,0,1); 
pNub = cdf('norm',+1,0,1);
p1sigma = pNub - pNlb

% Fdistribution
ndof1=9
ndof2=9
x=0.001:0.05:5;
Fc=cdf('F',x,ndof1,ndof2);
Fp=pdf('F',x,ndof1,ndof2);

% get bounds of F-statistic
Flb = icdf('F',pNlb,ndof1,ndof2)
Fub = icdf('F',pNub,ndof1,ndof2)


figure; 
subplot(2,1,1);
hold on;
plot(z,Nc,'r-');
plot(z,Np,'b-');
plot([-1 -1],[0, pNlb],'g-');
plot([+1 +1],[0, pNub],'g-');
xlabel('Z');
ylabel('P');
title('Standard normal (0, 1)');
legend('CDF','PDF','68%');

subplot(2,1,2);
hold on;
plot(x,Fc,'r-');
plot(x,Fp,'b-');
plot([Flb, Flb],[0, pNlb],'g-');
plot([Fub, Fub],[0, pNub],'g-');
xlabel('x');
ylabel('P');
title(sprintf('F (%d, %d)',ndof1, ndof2));
legend('CDF','PDF','68%');

