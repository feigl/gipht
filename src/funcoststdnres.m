function cost1 = funcoststdnres(DST,PST,TST)
%function cost1 = funcoststdnres(DST,PST,TST)
% 20160524 return chi2 statistic

narginchk(3, 3);

% field of residuals in radians
resids = funcostsdiff(DST,PST,TST);

% cost is sample standard deviation of weighted residuals
cost1=std(resids ./ DST.phasig); 

return;

