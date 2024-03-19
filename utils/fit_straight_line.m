function [pest, psig, tfit, ymod, ymodl, ymodu, mse] = fit_straight_line(time,yobs,ysig)
%function [pest, psig, tfit, ymod, ymodl, ymodu, mse] = fit_straight_line(time,yobs,ysig)
% fit a straight line to data using weighted least squares
% inputs:
%     time [arbitrary units]
%     yobs [arbitrary units]
%     ysig [arbitrary units]
% output:
%     pest: estimated values of parameters coeficients (intercept and slope
%     psig: corresponding uncertainties 
%     mse: mean squared error
% 20141125 Kurt Feigl
% 20220906 Kurt Feigl add help
% 20240319 Kurt Feigl square sigmas to obtain variances

% sanity check
n = numel(time);

% if ysig is not specified, then default to unweighted least squares
if nargin == 2
    ysig = ones(size(yobs));
end

if numel(yobs) ~= n || numel(ysig) ~= n
    error('miscount')
end

% prune out missing data
iok = 1:n;
iok = intersect(iok, find(isfinite(time)));
iok = intersect(iok, find(isfinite(yobs)));
iok = intersect(iok, find(isfinite(ysig)));

tfit = time(iok);
yobs = yobs(iok);
ysig = ysig(iok);
n = numel(iok);

tmid = mean(tfit);

% build design matrix
m = 2; % number of parameters
G = zeros(n,m);
for i=1:n
    j = 1; G(i,j) = 1.0;
    j = 2; G(i,j) = tfit(i) - tmid;
end

% solve using least squares
%[pest, psig, mse] = lscov(G, yobs, diag(ysig));
% 20240319 Kurt Feigl square sigmas to obtain variances
[pest, psig, mse] = lscov(G, yobs, diag(ysig .^2 ));

% calculate modeled values
ymod = G * pest;

% lower bounding envelope
ymodl = G * (pest - psig);

% upper bounding envelope
ymodu = G * (pest + psig);
       
return
end

