function bp=bperp(alon,alat,aheight,orbfile1,orbfile2)
% calculate perpendicular component of orbital separation ("baseline") vector between 
% two orbits, each at their nearest time of passage
% %
% usage: bp=bperp(alon,alat,aheight,orbfile11)
%
% alon    == geodetic longitude in degrees, eastward  reckoned positive
% alat    == geodetic latitude  in degrees, northward reckoned positive
% aheight == height above WGS84 ellipsoid in meters
% orbfile1 == ORB file in DIAPASON format, e.g., '12345.orb'.
%            see readorb.m for more details
% 
% Kurt Feigl 2011-JUN-27 test on Landers


if nargin ~= 5
   error  'Wrong number of arguments.'
end

if abs(alat) > 90
   error(sprintf('Impossible latitude: %f\n,',alat));
end

if alon > 360 | alon < -180
   error(sprintf('Impossible longitude: %f\n,',alon));
end

if aheight > 9000e3 | aheight < -100
   error(sprintf('Impossible height: %f\n,',aheight));
end

if fexist(orbfile1) ~= 1
   error(sprintf('Cannot find orbit file: %s\n,',orbfile1));
end

if fexist(orbfile2) ~= 1
   error(sprintf('Cannot find orbit file: %s\n,',orbfile2));
end

% position of point on ground
[xa, ya, za] = wgs84toxyz(alon, alat, aheight);

% state vectors of satellites
[xs1, ys1, zs1, xdot1, ydot1, zdot1, mjd1, sec1, orbnum1] = readorb(orbfile1);
[xs2, ys2, zs2, xdot2, ydot2, zdot2, mjd2, sec2, orbnum2] = readorb(orbfile2);

% position of satellite at nearest passage
[xn1, yn1, zn1, mjdn1, secn1] = nearestpassage(xs1, ys1, zs1, xa, ya, za, mjd1, sec1, xdot1, ydot1, zdot1);
[xn2, yn2, zn2, mjdn2, secn2] = nearestpassage(xs2, ys2, zs2, xa, ya, za, mjd2, sec2, xdot2, ydot2, zdot2);

% slave(2) w.r.t. master(1)
dx = xn2 - xn1;
dy = yn2 - yn1;
dz = zn2 - zn1;
dl = sqrt(dx^2 + dy^2 + dz^2)

% unit look vector from satellite to ground at master's nearest passage
ux = xa - xn1;
uy = xa - yn1;
uz = za - zn1;
ul = sqrt(ux^2 + uy^2 + uz^2);
ux = ux/ul;
uy = uy/ul;
uz = uz/ul;

% component parallel to unit vector = change in range 
bpara = [dx dy dz] * [ux uy uz]'
% component perpendicular to unit vector
bp    = NaN;

return








return
end
