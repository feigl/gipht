function d0 = circular_mean_deviation(theta, imode)
% circular_mean_deviation (theta, imode)
% given angles theta in radians,
% return circular mean deviation in radians
% imode == 0 do NOT remove mean
% imode == 1 REMOVE mean [default]
% from Mardia 1972
% 
% Kurt Feigl 2008-FEB-11
%
if nargin < 2
   imode = 1;
end
if imode == 1
   alpha = mean_direction(theta);
else
   alpha = 0.0;
end
% if numel(isfinite(theta)) > 1 && numel(isreal(theta)) > 1 ...
%         isfinite(alpha) == 1 && isreal(alpha) == 1
% 2012-JAN-10 isreal returns a scalar!!
if numel(isfinite(theta)) > 1 && isfinite(alpha) == 1
    if ~isreal(theta)
        warning('Detected imaginary values in theta. Taking real part.');
        theta = real(theta);
    end
    if ~isreal(alpha)
        warning('Detected imaginary values in alpha. Taking real part.');
        theta = real(alpha);
    end    
    theta_prime = mod(theta-alpha,2*pi);
    eta = min(theta_prime,2*pi-theta_prime);
    d0 = sum(eta)/numel(eta);
else
    d0 = NaN;
end
return
