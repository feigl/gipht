function l=loglikelihood3(p,fitfun,DST,PST,TST)
%function l=loglikelihood2(params,xdat,yobs,ysig,fitfun)
% modified from loglikelihood.m
% Example 11.4
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
%
% 20140113 Kurt Feigl
%

% observed values
yobs = DST.phaobs;

% set trial value of parameters
PST.p1 = p;

% calculate modeled values by evaluating fitting function
%ymod = feval(fitfun,PST,DST,TST);
ymod = feval(fitfun,DST,PST,TST);

% data uncertainty
ysig = DST.phasig;

% Compute the normalized residuals.
% size(yobs)
% size(ymod)
% size(ysig)
yresn=(yobs-ymod)./ysig;

% log likelihood 
l=(-1/2)*sum(yresn.^2);

return
end
