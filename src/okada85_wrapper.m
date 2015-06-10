function uENZ = okada85_wrapper(pg,xyobs,nu)
% wrapper around Okada85 routines
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% % Assigns input arguments
% e = varargin{1};
% n = varargin{2};
% depth = varargin{3};
% strike = varargin{4}*pi/180;	% converting STRIKE in radian
% dip = varargin{5}*pi/180;	% converting DIP in radian ('delta' in Okada's equations)
% L = varargin{6};
% W = varargin{7};
% rake = varargin{8}*pi/180;	% converting RAKE in radian
% slip = varargin{9};
% U3 = varargin{10};

    % pg(1)   pg(2)    pg(3)  pg(4)  pg(5)  pg(6)   pg(7)   pg(8)    pg(9)   pg(10)
    % Length - Width - Depth - Dip - Strike - East - North - Sslip - Dslip - Op
    % (m)     (m)    (m)    (deg) (deg)    (m)   (m)    (m)     (m)     (m)
    %                  upper   neg-
    %                  edge    ative

%	   E,N    : coordinates of observation points in a geographic referential 
%	            (East,North,Up) relative to fault centroid (units are described below)
%	   DEPTH  : depth of the fault centroid (DEPTH > 0)
%	   STRIKE : fault trace direction (0 to 360° relative to North), defined so 
%	            that the fault dips to the right side of the trace
%	   DIP    : angle between the fault and a horizontal plane (0 to 90°)
%	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
%	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
%	   RAKE   : direction the hanging wall moves during rupture, measured relative
%	            to the fault STRIKE (-180 to 180°).
%	   SLIP   : dislocation in RAKE direction (length unit)
%	   OPEN   : dislocation in tensile component (same unit as SLIP)

% 1 E#  32 Okada1_Length_in_m______________  2340.0000  2340.0000     0.0000        NaN     NaN  2340.0000
% 2 E#  33 Okada1_Width_in_m_______________  2780.0000  2780.0000     0.0000        NaN     NaN  2780.0000
% 3 E#  34 Okada1_Centroid_Depth_in_m______  3000.0000  3000.0000     0.0000        NaN     NaN  3000.0000
% 4 E#  35 Okada1_Dip_in_deg_______________    50.0000    50.0000     0.0000        NaN     NaN    50.0000
% 5 E#  36 Okada1_Strike_CCW_from_N_in_deg_   282.0000   282.0000     0.0000        NaN     NaN   282.0000
% 6 E#  37 Okada1_Easting_in_m_____________   507200.0   507200.0        0.0        NaN     NaN   507200.0
% 7 E#  38 Okada1_Northing_in_m____________  3802200.0  3802200.0        0.0        NaN     NaN  3802200.0
% 8 E#  39 Okada1_RL_Strike_Slip_in_m______  -1.99e-02  -1.99e-02   0.00e+00        NaN     NaN  -1.99e-02
% 9 E#  40 Okada1_Downdip_Slip_in_m________    -0.5630    -0.5630     0.0000        NaN     NaN    -0.5630
% 10 F#  41 Okada1_Tensile_Opening_in_m_____   0.00e+00   0.00e+00   0.00e+00        NaN     NaN   0.00e+00

% for i=1:numel(pg)
%     fprintf(1,'PG(%3d) = %10.2E\n',i,pg(i));
% end
STRIKE =    pg(5); % strike in degrees
DIP    =    pg(4); % dip in degrees
LENGTH =    pg(1);
WIDTH  =    pg(2);

if isfinite(pg(8)) == 1 && isfinite(pg(9)) == 1
    % in terms of downdip and right-lateral strike slip
%     RAKE   =    atan2(pg(8),pg(8)) * 180./pi; % rake in degrees
%     SLIP   =    hypot(pg(8),pg(9));           % slip magnitude in meters
    % parameterized in terms of rake
    RAKE = pg(9);
    SLIP = pg(8);
end
OPEN   =    pg(10);                       % tensile opening in meters

E =     xyobs(1,:) - pg(6);               % relative position of obs point wrt to fault centroid
N =     xyobs(2,:) - pg(7);               % relative position of obs point wrt to fault centroid
DEPTH = pg(3);                            % relative position of obs point wrt to fault centroid

NU    = nu;                               % Poisson's ratio

 
%
%	returns the following variables (same matrix size as E and N):
%	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
%	   uZE,uZN         : tilts (in rad * FACTOR)
%	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
 


%[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
[uE,uN,uZ] = okada85(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU);

% disp 'uE'; size(uE)

uENZ(1,:) = uE;
uENZ(2,:) = uN;
uENZ(3,:) = uZ;

return



