function cost1 = funcostrarc4(p1,DST,PST,TST)
%function cost1 = funcostrarc(DST,PST,TST)
% cost function for phase model is mean of dev = arc(obs,mod)
% 20160524 reduce number of arguments to 4

narginchk(4, 4);

PST.p1 = colvec(p1);

% field of costs (angular devation) in radians
costs = funcostsrarc(DST,PST,TST);

% 20140807 prune out NaN
iok = isfinite(costs);
costs = costs(iok);

% number of elements
n=numel(costs);

% average cost is L1 norm in radians
cost1=sum(colvec(costs))/n; 

% convert to cycles
cost1 = cost1/2.0/pi;

return;

