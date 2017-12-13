function cost1 = funcostrms(DST,PST,TST)
%function cost1 = funcostrms(DST,PST,TST)
% cost function for is sum of squared residuals 

% field of residuals in same units as data
resids = funcostsdiff(DST,PST,TST);

% cost is RMS of finite-valued residuals
cost1 = rms(resids(isfinite(resids)==1));

return;

