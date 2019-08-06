function cost = funcostrarcminlength(DST,PST,TST)
%function cost1 = funcostrarcminlength(DST,PST,TST)
% cost function for phase model is mean of dev = arc(obs,mod)
% 20160524 reduce number of arguments to 3
% 20190405 add a penalty for L1 roughness

narginchk(3, 3);

if isfield(PST,'beta') == 1
    beta = PST.beta;
else
    beta = 1.0;
end

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

% roughness as L1 norm of model (unit) vector
modelvec = PST.p1;
cost2 = sum(abs(diff(modelvec)))/norm(modelvec);

% total cost should be between 0 and 1
cost = (cost1 + beta*cost2)/(1.0 + beta); 

return;
end


