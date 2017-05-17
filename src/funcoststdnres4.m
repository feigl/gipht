function cost1 = funcoststdnres4(p1,DST,PST,TST)
%function cost1 = funcoststdnres(DST,PST,TST)
% 20160524 return chi2 statistic
% 20170516 take 4 arguments, first of which is parameter vector

narginchk(4, 4);

PST.p1 = colvec(p1);

% field of residuals 
resids = funcostsdiff(DST,PST,TST);

% cost is sample standard deviation of weighted residuals
cost1=std(resids ./ DST.phasig);

return;

