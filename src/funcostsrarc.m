function costs = funcostsrarc(fitfun,DST,PST,TST)
%function costs = funcostsrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials)
%function costs = funcostsrarc(p,fitfun,varargin)
%function costs = funcostsrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1)
% cost  function for phase model
%   p   == parameter
%   xyzm == easting, northing, up in m [3 rows] 
%   tepochs = unique times in decimal years [1 row]
%   bpest = orbital separation w.r.t. virutal orbit in meters [1 row]
%   dops = Doppler shifts w.r.t. virtual Dopple in PRF [1 row
%   DD = differencing matrix of 1s and 0s
%   unitv = unit vector
%   xd == observed (wrapped) phase gradient in radians [-pi,+pi]
%         on 1-byte signed integer
%   yd == dummy, for historical reasons
% return values in DN [0 to 127]
% for use with ANNEAL

% nargchk(10, 10, 10);
%
% if (size(xyzm) ~= size(xd))
%     error 'arguments wrong dimension'
% else

% evaluate fitting function at current value of parameters
% unwrapped phase in radians
%[uwm, dummyforpartials] = feval(fitfun,p,xyzm,tepochs,bpest,dops,DD,unitv,ippix1,ifast,partials);
%[uwm, dummyforpartials] = feval(fitfun,DST,PST,TST);
uwm = feval(fitfun,DST,PST,TST);

%     if (size(uwm) ~= size(xd))
%         size(uwm)
%         size(xd)
%         error 'size problem'
%     end

% phase gradient in radians on double
wrm=rwrapm(uwm);% call wrap function in Matlab

% arc function in radians on double
%costs = rarcm(wrm,xd);  % call matlab function
costs = rarcm(wrm,DST.phaobs);  % call matlab function

return;

