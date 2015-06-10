function [ha, bperp]=hamb(alon,alat,ahit,orbfilem,orbfiles)
% calculate altitude of ambiguity and perpendicular baseline from two
% orbital files
% 
% alon    == geodetic longitude in degrees, eastward  reckoned positive
% alat    == geodetic latitude  in degrees, northward reckoned positive
% ahit == height above WGS84 ellipsoid in meters
% orbfilem == master ORB file in DIAPASON format, e.g., '12345.orb'.
%            see readorb.m for more details
% orbfiles == slave ORB file in DIAPASON format, e.g., '54321.orb'.
%            see readorb.m for more details
% 
% 
% Kurt Feigl 2012-NOV-08
% 


if nargin ~= 4
   error  'Wrong number of arguments.'
end

if abs(alat) > 90
   error(sprintf('Impossible latitude: %f\n,',alat));
end

if alon > 360 || alon < -180
   error(sprintf('Impossible longitude: %f\n,',alon));
end

if ahit > 9000e3 || ahit < -100
   error(sprintf('Impossible height: %f\n,',ahit));
end

[xms, yms, zms, xmdot, ymdot, zmdot, mjdm, secm, orbnumm] = readorb(orbfilem);
[xss, yss, zss, xsdot, ysdot, zsdot, mjds, secs, orbnums] = readorb(orbfiles);

[xa, ya, za] = wgs84toxyz(alon, alat, ahit);

[slon, slat, srad] = carsph(xa, ya, za);

[xnm, ynm, znm, mjdnm, secnm, distnm, velonm] = nearestpassage(xms, yms, zms...
    , xa, ya, za...
    , mjdm, secm ...
    , xmdot, ymdot, zmdot);

[xns, yns, zns, mjdn, secns, distns, velons] = nearestpassage(xss, yss, zss...
    , xa, ya, za...
    , mjds, secs ...
    , xsdot, ysdot, zsdot);


% from ground to satellite
dx = xnm - xa;
dy = ynm - ya;
dz = znm - za;

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
