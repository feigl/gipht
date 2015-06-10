function [U, D, N, R, H, A, V, mjdn, secn]=orbvectors2(tlon,tlat,ahit ...
    ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum)
% function [U, D, N, R, H, A, V, mjdn, secn]=orbvectors2(tlon,tlat,ahit ...
%   ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum)
% calculate orbital quantities
%
% inputs
%    tlon    == target geodetic longitude in degrees, eastward  reckoned positive
%    tlat    == target geodetic latitude  in degrees, northward reckoned positive
%    ahit    == target height above WGS84 ellipsoid in meters
%    orbfile == ORB file in DIAPASON format, e.g., '12345.orb'.
%            see readorb.m for more details
% OUTPUTS
%
%  mjdn             == epoch of satellite state vector at nearest passage in Modified Julian Day
%  secn             == epoch of satellite state vector at nearest passage in seconds of day
%
% U = [ue, un, uu] == unit vector pointing from ground to satellite 
%    in local [Eastward, Northward, Upward] 
%    such that range change = -1 * [disE disN disU] * [ue un uu]^T  
%    Be sure to remember the -1 in the equation above!
%
% D = [dx, dy, dz] == Differential Position (in meters)
%                     of Target w.r.t. Satellite at Nearest Passage
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
% 
% N = [nx, ny, nz] == Position (in meters) of Satellite at Nearest Passage
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
%
% R = [rx, ry, rz] == Velocity in (m/s) of Satellite at Nearest Passage
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
%
% H = [hx, hy, hz] == horizontal component of look vector 
%                     along line of sight from satellite to target
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
%
% A = [ax, ay, az] == along-track component of look vector 
%                     along line of sight from satellite to target
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
%
% V = [vx, vy, vz] == vertical component of look vector 
%                     along line of sight from satellite to target
%                     in Cartesian, Geocentric, Inertial Frame [X, Y, Z]
%
% Kurt Feigl
% 2011-JUL-01

if nargin ~= 12
   error  'Wrong number of arguments.'
end

if abs(tlat) > 90
   error(sprintf('Impossible latitude: %f\n,',tlat));
end

if tlon > 360 || tlon < -180
   error(sprintf('Impossible longitude: %f\n,',tlon));
end

if ahit > 9000e3 || ahit < -100
   error(sprintf('Impossible height: %f\n,',ahit));
end

% position of target
[xt, yt, zt] = wgs84toxyz(tlon, tlat, ahit);
T = [xt; yt; zt]; % position of target written as column vector

% position and velocity of satellite
%[xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum] = readorb(orbfile);
 
% position and velocity of satellite at nearest passage
%tstart1 = tic; 
[xn, yn, zn, un, vn, wn, mjdn, secn] = nearestpassage (xs, ys, zs, xt, yt, zt, mjd, secs, xdot, ydot, zdot);
%telapsed1 = toc(tstart1)
%tstart2 = tic; 
%[xn, yn, zn, un, vn, wn, mjdn, secn] = nearestpassage2(xs, ys, zs, xt, yt, zt, mjd, secs, xdot, ydot, zdot);
%telapsed2 = toc(tstart2)

N = [xn; yn; zn]; % position of satellite at nearest passage written as column vector
R = [un; vn; wn]; % velocity of satellite at nearest passage written as column vector

% differential position of satellite w.r.t. target 
D = N - T;
dn = sqrt(abs(D(1))^2 + abs(D(2))^2 + abs(D(3))^2); % avoid overflow
D = D/dn;

% rotate to local coordinates tangent to ellipsoid at target
[ue, un, uu] = xyz2local2 (tlon, tlat, ahit, D(1), D(2), D(3));
% normalize to unit length
ud = sqrt(abs(ue)^2 + abs(un)^2 + abs(uu)^2); % avoid overflow
U = [ue; un; uu]/ud; % write as a column vector

% Along-track component is parallel to velocity
A = R / sqrt(abs(R(1))^2 + abs(R(2))^2 + abs(R(3))^2); % avoid overflow

% "Vertical" component is parallel to radius
V = N / sqrt(abs(N(1))^2 + abs(N(2))^2 + abs(N(3))^2); % avoid overflow

% Horizontal component is across track
%H = cross(V,A);
H = cross(A,V); % use same sign convention as Diapason alt_ambiguity
H = H / sqrt(abs(H(1))^2 + abs(H(2))^2 + abs(H(3))^2); % avoid overflow

return
end
