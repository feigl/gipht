function [de, dn, du] = xyz2local2 (alon, alat, ahit, x, y, z)
% rotate x, y, z to e,n,u 
% tangent  at ellipsoidal (WGS84) coordinates 
%[alon(in degrees), alat(in degrees), ahit(in meters)]

[xa, ya, za] = wgs84toxyz(alon, alat, ahit);

[geod_pos, R] = xyz_to_geod([xa, ya, za]);

% disp '              latitude  in degrees'; alat1 = 180.0*geod_pos(1)/pi
% disp '              longitude in degrees'; alon1 = 180.0*geod_pos(2)/pi
% disp '              height in meters    '; ahit1 =       geod_pos(3)
% 
% disp 'rotation matrix R such that R * [colat long height] = [X Y Z]'
% R
% 

DNEU =  R * [x y z]';

de =       DNEU(2); % east
dn =       DNEU(1); % latitude
du =       DNEU(3); % upward


return;
