function cost1 = funcostrmsnres(DST,PST,TST)
% return objective function as RMS of weighted residuals
% 20171213
%narginchk(3, 3);

% field of residuals in same units as data
resids = funcostsdiff(DST,PST,TST);

% cost is root mean square of weighted residuals
cost1=rms(resids ./ DST.phasig);

return;

