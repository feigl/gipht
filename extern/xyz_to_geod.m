function [geod_pos, rot_mat] = xyz_to_geod(xyz_pos)
% adapted from GAMIT ROUTINE of the same name - 2011-06-28 Kurt Feigl
%
% *     Routine to compute the geodetic latitudes, longitude and
% *     elipsoidal height from the cartesian XYZ coordinates of
% *     the site.  The algorithm used is adopted from:
% *     Heiskanen W A and H. Moritz, Physical Geodesy, W.H.Freeman
% *     and company, San Francisco, 364, 1967. Chapter 5.3.
% *
% *     The routine also returns the tranformation between XYZ and
% *     local NEU system.
% *
% *     RESTRICTION: User must supply the XYZ coordinates.  If these
% *     are not available the can be computed from geod_pos using:
% *
% *     X = (N+h) cos(phi) cos(lambda)
% *     Y = (N+h) cos(phi) sin(lambda)
% *     Z = [(1-e^2)N+h]  sin(phi)
% *
% *     where  e^2 = 2f-f^2; f is flattening
% *            N^2 = a^/[1-(e*sin(phi))^2]
% *     N is East-West Radius of curvature.
% *
% *   dlatdN      - Derivative of latitudes with respect to the
% *               - (horizontal) North component of position
% *               - This derivative should not be confused with
% *               - dNdlat (same symbol N but different meanings)
% *   dNdlat      - the derivative of the curvature (rad_curve)
% *               - with respect to latitudes.
% *   dXdlat      - Derivative of the X coordinate wrt latitudes
% *   dYdlat      - Derivative of the Y coordinate wrt latitudes
% *   dZdlat      - Derivative of the Z coordinate wrt latitudes
%
% *   eccsq       - eccentricity squared computed from flattening
% *   equ_rad     - radial distance of site from rotation axis.
% *   geod_pos(3) - the geodetic coordinates.  The order of the
% *               - values is:
% *   XXXXXXXXXX  - (1) - Geodetic co-latitudes (rads)
% *   MOD         - (1) - Geodetic   latitudes (rads)
% *               - (2) - Geodetic longitude (positive east) (rad)
% *               - (3) - ellipsoidal height (m)
% *   lat_i       - approximate latitudes used in iteration
% *   lat_p       - approximate latitudes from previous iteration
% *   long        - longitude of the site.
% *   h_i         - approximate height used in iteration
% *   h_p         - approximate height from previous iteration
% *   rad_lat     - radius to be used in computing the geodetic
% *               - latitudes
% *   rad_curve    - radius of curvature at the site (N above)
% *   rot_mat(3,3)    - Transformation matrix XYZ and NEU where
% *               - N is -colatitudes, E is longitude and U is
% *               - along ellipsoidal height
% *   xyz_pos(3) - XYZ coordinates of the site (m)
%
% *   tolerance   - toleranace for convergence on geodetic
% *               - coordinates.
%
%       real*8 dlatdN, dNdlat, dXdlat, dYdlat, dZdlat, eccsq, equ_rad,
%      .    geod_pos(3), lat_i, lat_p, long, h_i, h_p, rad_lat,
%      .    rad_curve, rot_mat(3,3), xyz_pos(3), tolerance,
%      .    earth_rad, earth_flat, pi
%
%       integer*4 len,rcpar
%
%       character*80 prog_name
%
%  WGS-84 parameters for the elliposoid (a and 1/f)
earth_rad     = 6378137.D0;
earth_flat    = 1.d0/298.257222101D0;


%
% *   niter       - Number of iterations.  EXITS if more than 50
% *                 iterations are required.  (Removes problem if
% *                 coordinates are NaN in IEEE floating point)
%
%       integer*4 niter
%
% *   converged   - Indicate values have converged.
%
%       logical converged
%
%       data
% Converge to 0.1 mm
tolerance  = 0.0001d0;

% Start, get the latitudes and height by iteration.
equ_rad = sqrt( xyz_pos(1)^2 + xyz_pos(2)^2 );
eccsq   = 2.d0*earth_flat - earth_flat^2;

%  Set previous iteration values
lat_p = atan2( xyz_pos(3), equ_rad);
h_p   = 0.d0;

converged = 0;
niter = 0;

geod_pos = nan(size(xyz_pos));
while (converged == 0)    
    %   Get radius of curvature using previous latitudes estimate
    rad_curve = earth_rad / sqrt(1.d0 - eccsq*(sin(lat_p)^2));
    rad_lat  = equ_rad    * (1.d0 - eccsq*rad_curve/(rad_curve+h_p));
    
    lat_i = atan2( xyz_pos(3), rad_lat);
    if( abs(lat_i) < pi/4.0 )
        %   Use cos lat formula
        h_i   = equ_rad/cos(lat_i) - rad_curve;        
    else
        %  Use sin lat formula
        h_i   = xyz_pos(3)/sin(lat_i) - (1.d0-eccsq)*rad_curve;
    end
       
    %         Check for convergence
    if( abs(h_i-h_p) < tolerance && abs(lat_i-lat_p)*rad_curve < tolerance )
        converged = 1;
    end
    
    %         Check for two many iterations
    niter = niter + 1;
    if( niter > 50 ) 
        warning('FATAL lib/xyz_to_geod Convergence failure in computing geodetic coordinates');
        return;
    end
    
    % Save the latest values
    h_p   = h_i;
    lat_p = lat_i;    
end % iterating for latitudes and height

% Save the final values
long = atan2(xyz_pos(2), xyz_pos(1));

% colatitudes
%geod_pos(1) = pi/2.d0 - lat_i;
% Mod: 2011-JUN-28 Kurt Feigl
% latitudes
geod_pos(1) = lat_i;


% Add 2*pi
if( long <0.0 ) 
    geod_pos(2) = 2.0*pi + long;
else
    geod_pos(2) = long;
end

geod_pos(3) = h_i;

%Now do the transformation between XYZ and NEU

%     Compute derivative of the radius of curvature with to latitudes
dNdlat = 2*earth_rad*eccsq *sin(lat_i)*cos(lat_i) / ...
    sqrt( (1.d0-eccsq* sin(lat_i)^2)^3 );

%     Now do NORTH component  (First do derivate wrt latitudes, then compute
%     the Northing derivative)
dXdlat  = (dNdlat * cos(lat_i) - ...
    (rad_curve+h_i)* sin(lat_i) )*cos(long);
dYdlat  = (dNdlat * cos(lat_i) - ...
    (rad_curve+h_i)* sin(lat_i) )*sin(long);
dZdlat  = (1.d0-eccsq)*dNdlat         *sin(lat_i) + ...
    ((1.d0-eccsq)*rad_curve + h_i)*cos(lat_i);

dlatdN  = sqrt(dXdlat^2 + dYdlat^2 + dZdlat^2 );

%     Now do rotation matrix
rot_mat(1,1) = dXdlat / dlatdN;  % dX/dN
rot_mat(1,2) = dYdlat / dlatdN;  % dY/dN
rot_mat(1,3) = dZdlat / dlatdN;  % dZ/dN

%     Now do EAST component
rot_mat(2,1) = -sin(long);
rot_mat(2,2) =  cos(long);
rot_mat(2,3) =  0.d0;

%     Now do UP component
rot_mat(3,1) = cos(lat_i) * cos(long); 
rot_mat(3,2) = cos(lat_i) * sin(long);
rot_mat(3,3) = sin(lat_i);

% Thats all
return
end

