function [xn, yn, zn, un, vn, wn, mjdn, secn] = nearestpassage(xs, ys, zs, xt, yt, zt, mjd, secs, xdot, ydot, zdot)
%function [xn, yn, zn, un, vn, wn, mjdn, secn] = nearestpassage(xs, ys, zs,
%xt, yt, zt, mjd, secs, xdot, ydot, zdot)
% find epoch nearest to target at xt, yt, zt
%
% INPUTS:
%
%  mjd              == epochs of satellite state vectors in Modified Julian Day
%  secs              == epochs of satellite state vectors in seconds of day
%  xt, yt, zt       == Position Vectors of Target
%                      in Cartesian Geocentric X, Y, Z meters
%  xs, ys, zs       == Position Vectors of Satellite 
%                      in Cartesian Geocentric X, Y, Z meters
%  xdot, ydot, zdot == Velocity Vectors of Satellite 
%                      in Cartesian Geocentric X, Y, Z meters/second
%  
% OUTPUTS:
%
%  mjdn             == epoch of satellite state vector at nearest passage in Modified Julian Day
%  secn             == epoch of satellite state vector at nearest passage in seconds of day
%  xn, yn, zn       == Position of Satellite at nearest passage
%                      in Cartesian Geocentric X, Y, Z meters
%  un, vn, wn,      == Velocity Vector of Satellite at nearest passage
%                      in Cartesian Geocentric X, Y, Z meters/second
%  Kurt Feigl 2011-JUL-02

mjd0 = min(mjd);

t = 86400.0 * (mjd - mjd0) + secs;

n = numel(xs);

if numel(ys) ~= n || numel(zs) ~= n
    error('dimension problem');
end

if exist('xdot','var') ~= 1
    %xdot = zeros(size(xs));
    xdot = diff(xs) ./ diff(t);
    xdot(end+1) = NaN;
end
if exist('ydot','var') ~= 1
    %ydot = zeros(size(ys));
    ydot = diff(ys) ./ diff(t);
    ydot(end+1) = NaN;
end
if exist('zdot','var') ~= 1
    %zdot = zeros(size(zs));
    zdot = diff(zs) ./ diff(t);
    zdot(end+1) = NaN;
end

% interpolate the state vectors
%dt = 0.001; % step size in seconds

i1 = 1;
i2 = n;
tstep = mean(diff(t));
dt = tstep/10.0;
while dt > 0.001
    % distance between satellite at [xs, ys, zs] and target at [xt, yt, zt]
     dist = sqrt ( (xs - xt).^2 + (ys - yt).^2 + (zs - zt).^2 );
    
%     fprintf (1,'i, t, dist, dt = %12.4f\n',dt);
%     for i=1:n
%       fprintf (1,'%6d %12.4f %12.4f\n',i,t(i),dist(i));
%     end
    
    % pointer to minimum distance
    imin = find(abs(dist-min(dist)) < 1.0e-6);
    if numel(imin) > 0
        imin = imin(1);
    else
        break;
    end
    
    % recover quantities at time of nearest passage
    xn = xs(imin);
    yn = ys(imin);
    zn = zs(imin);
    un = xdot(imin);
    vn = ydot(imin);
    wn = zdot(imin);
    distn = dist(imin);
    secn = mod(t(imin),86400.0);
%    mjdn = floor(t(imin)/86400.0) + min(mjd);
    mjdn = floor(t(imin)/86400.0) + mjd0;
    
%     % draw a picture
%     figure; hold on;
%     plot(t,dist,'k.-');
%     plot(t(imin),dist(imin),'ro');
%     plot([min(t)   max(t)],[dist(imin) dist(imin)],'r-');
%     plot([t(imin) t(imin)],[ min(dist)  max(dist)],'r-');
%     xlabel('time (s)');
%     ylabel('distance in meters');
%     title(sprintf('Nearest passage is %#20.4f meters at MJD %16d and %#20.4f seconds\n',dist(imin),mjdn,secn));
    
    % % Interpolate using Matlab splines
    i1 = max([imin-2,1]);
    i2 = min([imin+2,n]);

    ti    = [       t(i1):dt:t(i2)];
    xi    = interp1(t(i1:i2),xs(i1:i2)  ,ti,'spline');
    yi    = interp1(t(i1:i2),ys(i1:i2)  ,ti,'spline');
    zi    = interp1(t(i1:i2),zs(i1:i2)  ,ti,'spline');
    xdoti = interp1(t(i1:i2),xdot(i1:i2),ti,'spline');
    ydoti = interp1(t(i1:i2),ydot(i1:i2),ti,'spline');
    zdoti = interp1(t(i1:i2),zdot(i1:i2),ti,'spline');
%     mjdi  = interp1(t(i1:i2),mjd(i1:i2),ti,'spline');
    
    % set up for next iteration
    xs = xi;
    ys = yi;
    zs = zi;
    xdot = xdoti;
    ydot = ydoti;
    zdot = zdoti;
    t    = ti;
    n = numel(xs);
    dt = dt/10;
%     mjd = mjdi;
end
return;

