function l=loglikelihood2(params,xdat,yobs,ysig,fitfun)
%function l=loglikelihood2(params,xdat,yobs,ysig,fitfun)
% modified from loglikelihood.m
% Example 11.4
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
%
% 20140107 Kurt Feigl
%

% modeled values
ymod = feval(fitfun,params,xdat);

% Compute the normalized residuals.
% size(yobs)
% size(ymod)
% size(ysig)
yresn=(yobs-ymod)./ysig;

% log likelihood 
l=(-1/2)*sum(yresn.^2);

return
end
