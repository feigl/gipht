function [slon, slat, srad] = carsph(x,y,z)
% convert cartesian X,Y,Z in meters to spherical coordinates in degrees
srad = sqrt(x.^2 + y.^2 +z.^2);
slon = 180.0*(atan2(y,x))/pi;
slat = 180.0*(asin(z./srad))/pi;
return;
