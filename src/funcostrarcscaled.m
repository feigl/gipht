function cost1 = funcostrarcscaled(DST,PST,TST)
%function cost1 = funcostrarcscaled(DST,PST,TST)
% 2014-JAN-08 angular deviations by measurement uncertainty
% 20160524 change number of arguments

narginchk(3, 3);


% field of costs (angular devation) in radians
costs = funcostsrarc(DST,PST,TST);

% normalize by uncertainty
costs = costs ./ DST.phasig;

% number of elements
n=numel(costs);

% average cost is L1 norm in radians
cost1=sum(colvec(costs))/n; 

% convert to weighted cycles
cost1 = cost1/2.0/pi;

return;

