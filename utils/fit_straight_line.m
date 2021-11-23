function [pest, psig, tfit, ymod, ymodl, ymodu, mse] = fit_straight_line(varargin)
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
% example
%    time = [1995:1:2020]';
%    yobs = 0.01 * time - 19.95 + 0.0005*time .* randn(size(time));
%    ysig = 2.0 * ones(size(time));
%    [pest, psig, tfit, ymod, ymodl, ymodu, mse] = fit_straight_line(time,yobs,ysig);
%    pest
%    psig
%    figure; hold on;
%    errorbar(time,yobs,ysig,'ro');
%    plot(tfit,ymod,'k-');
% 20150614 Kurt Feigl
% 

narginchk(2,3);

time = varargin{1};
yobs = varargin{2};

% sanity check
n = numel(time);

% if ysig is not specified, then default to unweighted least squares
if nargin == 3
    ysig = varargin{3};
else
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

% solve using weighted least squares
%[pest, psig, mse] = lscov(G, yobs, diag(ysig)); 
% use variance 20160606
[pest, psig, mse] = lscov(G, yobs, diag(ysig.^2));

% calculate modeled values
ymod = G * pest;

% lower bounding envelope
ymodl = G * (pest - psig);

% upper bounding envelope
ymodu = G * (pest + psig);
        
    
return

end

