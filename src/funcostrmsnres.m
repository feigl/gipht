function cost1 = funcostrmsnres(DST,PST,TST)
% return objective function as RMS of weighted residuals
% 20170416
%narginchk(3, 3);

% field of residuals in radians
resids = funcostsdiff(DST,PST,TST);

% cost is sample standard deviation of weighted residuals
cost1=rms(resids ./ DST.phasig);

return;

