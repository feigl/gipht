function costs = funcostsdiff(DST,PST,TST)
%function costs = funcostsdiff(DST,PST,TST)
% 20160524 return residual differences = observed minus calculated
narginchk(3, 3);

% evaluate fitting function 
fitfun = str2func(PST.fitfun);
uwm = feval(fitfun,DST,PST,TST);

% residual 
costs = DST.phaobs - uwm; 

iok = find(isfinite(costs)==1);

return;

