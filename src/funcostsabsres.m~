function costs = funcostsabsres(fitfun,DST,PST,TST)
%function costs = funcostsabsres(fitfun,DST,PST,TST)
% return absolute value of residuals
% evaluate fitting function at current value of parameters
% modeled value 
modval = feval(fitfun,DST,PST,TST);

costs = abs(modval-DST.phaobs);  

return;

