function cost1 = funcostrarc(DST,PST,TST)
%function cost1 = funcostrarc(DST,PST,TST)
% cost function for phase model is mean of dev = arc(obs,mod)
% 20160524 reduce number of arguments to 3

narginchk(3, 3);

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
end
