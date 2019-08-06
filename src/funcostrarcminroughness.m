function cost = funcostrarcminroughness(DST,PST,TST)
%function cost = funcostrarcminroughness(DST,PST,TST)
% cost function for phase model is mean of dev = arc(obs,mod)
% 20160524 reduce number of arguments to 3
% 20190405 add a penalty for L1 roughness

narginchk(3, 3);

if isfield(PST,'beta') == 1
    beta = PST.beta;
elseif isfield(PST,'Roughness_Smoothing_Beta_Dimless') == 1
    beta = PST.Roughness_Smoothing_Beta_Dimless; 
else
    %beta = 0.25; % will weight equally with random phase noise
    %beta = 0.125; % half as weighty 
    beta =  0.06125; % quarter as weighty 
    beta =  0.03; % less weighty 
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

% extract part of parameter vector containing slips
slipvecs = PST.p1(TST.koffset:end);
%numel(slipvecs)

% roughness as L1 norm of model (unit) vector
%cost2 = sum(abs(diff(slipvecs)))/norm(slipvecs);

% roughness as gradient
slipgrid = reshape(slipvecs,TST.faultnrows,TST.faultncols);
[gradslipx,gradslipy] = gradient(slipgrid);
cost2 = hypot(gradslipx,gradslipy);
cost2 = mean(colvec(cost2))/norm(colvec(cost2));

% total cost should be between 0 and 1
if isfinite(cost1) && isfinite(cost2)
    cost = (cost1 + beta*cost2)/(1.0 + beta);
elseif isfinite(cost1)
    cost = cost1;
else
    cost = nan;
end

return;
end


