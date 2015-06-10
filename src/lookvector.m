function [ue, un, uu, mjdn, secn, distn, velon, incid]=lookvector(alon,alat,ahit,orbfile)
% calculate unit look vector from ground to satellite 
% 
% such that range change = -1 * [disE disN disU] * [ue un uu]^T  
%
% Be sure to remember the -1 in the equation above
%
% usage: [ue, un, uu]=lookvector(alon,alat,ahit,orbfile)
%
% alon    == geodetic longitude in degrees, eastward  reckoned positive
% alat    == geodetic latitude  in degrees, northward reckoned positive
% ahit == height above WGS84 ellipsoid in meters
% orbfile == ORB file in DIAPASON format, e.g., '12345.orb'.
%            see readorb.m for more details
% [ue,un,uu] = lookvector(-116.91,34.36,100.,'/data/winsar/orbits/delft/4051.orb')
% ue =
%    0.372268867695445
% un =
%   -0.076427423094856
% uu =
%    0.924972831570653
%
% compare to results from alt_ambiguite in Diapason
%
% ue0 = 0.324203
% un0 =-0.070934
% uu0 = 0.943324
% 
% 
% Kurt Feigl 2007-APR-29
% 2009-JUL-22 check arguments
% 2011-JUN-27 test on Landers


if nargin ~= 4
   error  'Wrong number of arguments.'
end

if abs(alat) > 90
   error(sprintf('Impossible latitude: %f\n,',alat));
end

if alon > 360 | alon < -180
   error(sprintf('Impossible longitude: %f\n,',alon));
end

if ahit > 9000e3 | ahit < -100
   error(sprintf('Impossible height: %f\n,',ahit));
end

[xs, ys, zs, xdot, ydot, zdot, mjd, sec, orbnum] = readorb(orbfile);

[xa, ya, za] = wgs84toxyz(alon, alat, ahit);

[slon, slat, srad] = carsph(xa, ya, za);

[xn, yn, zn, mjdn, secn, distn, velon] = nearestpassage(xs, ys, zs...
    , xa, ya, za...
    , mjd, sec ...
    , xdot, ydot, zdot);

% from ground to satellite
dx = xn - xa;
dy = yn - ya;
dz = zn - za;

%dn = sqrt(dx^2 + dy^2 + dz^2);
%dn = norm([dx, dy, dz]);
dn = sqrt(abs(dx)^2 + abs(dy)^2 + abs(dz)^2);
dx = dx/dn;
dy = dy/dn;
dz = dz/dn;

% spherical
%[ue1, un1, uu1] = xyz2local (slon, slat, dx, dy, dz);

% ellipsoidal
[ue2, un2, uu2] = xyz2local2 (alon, alat, ahit, dx, dy, dz);

% compare result to Diapason
%disp 'differences in [E, N, U]'; [ue1-ue2, un1-un2, uu1-uu2]
% spherical
% Difference in dimensionless unit vector            
%   -0.000000750779595   0.002937017231314   0.000216837449437
% ue = ue1;
% un = un1;
% uu = uu1;
% ellipsoidal
% Difference in dimensionless unit vector            
%   -0.000000750779596   0.002944368594867   0.000000240259079
ue = ue2;
un = un2;
uu = uu2;

% normalize
ud = norm([ue, un, uu]);
ue = ue/ud;
un = un/ud;
uu = uu/ud;

% incidence angle from local vertical in radians
incid = acos(uu);

if incid < 0. || incid > pi/2.
    fprintf(1,'Unrealistic incidence angle %12.4g\n',incid);
    alon
    alat
    ahit
    warning('Unrealistic incidence angle');
end
return
end
