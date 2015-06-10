function [xa, ya, za] = wgs84toxyz(alon, alat, hght, finv, semi)
%function [xa, ya, za] = wgs84toxyz(alon, alat, hght, finv, semi)
% convert geodetic coordinates in degrees to cartesian coordinates in meters
% Roots: TFORM.F from GAMIT

% parameters for WGS84
if exist('finv','var') ~= 1
    %finv = 298.2572;
    finv = 298.257222101D0;
end
if exist('semi','var') ~= 1
    semi = 6378137.0;
end


% parameters for spherical earth
% finv = Inf;
% semi = 6.371525288789591e+06;

twopi= 2.0*pi;
f= 1.0/finv;
e2= 2.0*f - f*f;

sinlat= sin(twopi*alat/360.0);
coslat= cos(twopi*alat/360.0);
sinlon= sin(twopi*alon/360.0);
coslon= cos(twopi*alon/360.0);
curvn= semi./(sqrt(1.0-e2*sinlat.*sinlat));
 
xa= (curvn+hght).*coslat.*coslon;
ya= (curvn+hght).*coslat.*sinlon;
za= (curvn.*(1.0-e2)+hght).*sinlat;
 
return;
