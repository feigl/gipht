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
% 2023/01/24 update to handle datetimes

narginchk(2,3);

if isdatetime(varargin{1})
    handle_datetimes=true;
    t=varargin{1};
    fmt=t.Format;
    tz=t.TimeZone;
    time=datenum(t);
else
    time = varargin{1};
    handle_datetimes=false;
end

yobs = varargin{2};

% sanity check
n = numel(time);

if nargin == 3
     ysig = varargin{3};
end
% if ysig is not specified, then default to unweighted least squares
% 2023/01/09 - this is wrong because it sets sigma to 1 (unit)
% else
%     ysig = ones(size(yobs));
% end

if nargin == 3
    if numel(yobs) ~= n || numel(ysig) ~= n
        error('miscount')
    end
else
    if numel(yobs) ~= n
        error('miscount')
    end
end

% prune out missing data
iok = 1:n;
iok = intersect(iok, find(isfinite(time)));
iok = intersect(iok, find(isfinite(yobs)));
if nargin == 3
    iok = intersect(iok, find(isfinite(ysig)));
end
n = numel(iok);

tfit = time(iok);
yobs = yobs(iok);
if nargin == 3
    ysig = ysig(iok);
end


tmid = mean(tfit);

% build design matrix
m = 2; % number of parameters
G = zeros(n,m);
for i=1:n
    j = 1; G(i,j) = 1.0;
    j = 2; G(i,j) = tfit(i) - tmid;
end

cond(G)

% solve using weighted least squares
%[pest, psig, mse] = lscov(G, yobs, diag(ysig)); 
% use variance 20160606
% corrected 2023/01/09 Kurt Feigl
if nargin == 3
    [pest, psig, mse] = lscov(G, yobs, diag(ysig.^2));
else
    [pest, psig, mse] = lscov(G, yobs);
end

%    lscov assumes that the covariance matrix of B is known only up to a
%     scale factor.  MSE is an estimate of that unknown scale factor, and
%     lscov scales the outputs S and STDX appropriately.  However, if V is
%     known to be exactly the covariance matrix of B, then that scaling is
%     unnecessary.  To get the appropriate estimates in this case, you should
%     rescale S and STDX by 1/MSE and sqrt(1/MSE), respectively.
% 
% calculate modeled values
ymod = G * pest;

% lower bounding envelope
ymodl = G * (pest - psig);

% upper bounding envelope
ymodu = G * (pest + psig);

% convert tfit back to datetime
if handle_datetimes
    tfit=datetime(datestr(tfit));
    tfit.Format=fmt;
    tfit.TimeZone=tz;
end

% figure;hold on;
% plot(tfit,yobs,'r+')
% plot(tfit,ymodl,'bo-')

        
return

end

