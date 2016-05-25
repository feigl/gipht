function costs = funcostsdiff(DST,PST,TST)
%function costs = funcostsdiff(DST,PST,TST)
% 20160524 return residual differences = observed minus calculated
narginchk(3, 3);

% evaluate fitting function in radians
fitfun = str2func(PST.fitfun);
uwm = feval(fitfun,DST,PST,TST);

% residual in radians
costs = DST.phaobs - uwm;  

return;

