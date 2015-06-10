function cost1 = funcostrarc(p,fitfun,DST,PST,TST)
%function cost1 = funcostrarc(p,fitfun,DST,PST,TST)
% cost function for phase model is mean of dev = arc(obs,mod)

nargchk(5, 5, nargin);

if numel(p) ~= numel(PST.p0)
    error(sprintf('Dimension mismatch %d %d\n',numel(p),numel(PST.p0)));
end
PST.p1 = p;
% field of costs in DN [0,127]
%costs = funcostsiarc1(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1);
%costs = funcostsiarc1(p,fitfun,varargin{:});
%costs = funcostsiarc1(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
% field of costs on DN [0, 
%costs = funcostsrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
costs = funcostsrarc(fitfun,DST,PST,TST);
%size(costs)

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

