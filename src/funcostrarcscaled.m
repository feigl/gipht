function cost1 = funcostrarcscaled(p,fitfun,DST,PST,TST)
% cost  function for phase model
%   p   == parameter
% 2014-JAN-08 angular deviations by measurement uncertainty
%
% for use with ANNEAL

nargchk(5, 5, nargin);

if numel(p) ~= numel(PST.p0)
    error(sprintf('Dimension mismatch %d %d\n',numel(p),numel(PST.p0)));
end
PST.p1 = p;

% field of costs in radians
costs = funcostsrarc(fitfun,DST,PST,TST);

% normalize by uncertainty
costs = costs ./ DST.phasig;

% number of elements
n=numel(costs);

% average cost is L1 norm in radians
cost1=sum(colvec(costs))/n; 

% convert to weighted cycles
cost1 = cost1/2.0/pi;

return;

