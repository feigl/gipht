function cost1 = funcostrms(p,fitfun,DST,PST,TST)
%function cost1 = funcostrms(p,fitfun,DST,PST,TST)
% cost function for is sum of squared residuals 

nargchk(5, 5, nargin);

if numel(p) ~= numel(PST.p0)
    error(sprintf('Dimension mismatch %d %d\n',numel(p),numel(PST.p0)));
end
PST.p1 = p;

% field of absolute values of residuals
costs = funcostsabsres(fitfun,DST,PST,TST);

% 20140807 prune out NaN
iok = isfinite(costs);
costs = colvec(costs(iok));

% number of elements
n=numel(costs);

% cost is RMS L2 norm 
cost1=sqrt((costs'*costs)/n); 

% convert to cycles
%cost1 = cost1/2.0/pi;

return;

