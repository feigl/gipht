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

% F#  31 Okada1_Length_in_m______________  4308.0738  4308.0738     0.0000        NaN     NaN   424.3164
% F#  32 Okada1_Width_in_m_______________  5035.3854  5035.3854     0.0000        NaN     NaN   314.9414
% F#  33 Okada1_Centroid_Depth_in_m______  5999.9760  5999.9760     0.0000        NaN     NaN    82.5195
% F#  34 Okada1_Dip_in_deg_______________    65.9941    65.9941     0.0000        NaN     NaN     4.1699
% F#  35 Okada1_Strike_CW_from_N_in_deg__   206.0934   206.0934     0.0000        NaN     NaN     7.1875
% F#  36 Okada1_Centroid_Easting_in_m____   530996.9   530996.9        0.0        NaN     NaN      183.6
% F#  37 Okada1_Centroid_Northing_in_m___  4242682.7  4242682.7        0.0        NaN     NaN      431.2
% F#  38 Okada1_Coplanar_slip_in_m_______     0.2919     0.2919     0.0000        NaN     NaN     0.0447
% F#  39 Okada1_Rake_in_deg_CCW__________   -28.7569   -28.7569     0.0000        NaN     NaN     4.8975
% F#  40 Okada1_Tensile_Opening_in_m_____   0.00e+00   0.00e+00   0.00e+00        NaN     NaN   0.00e+00

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

if abs(OPEN) < 0.01 && abs(SLIP) < 0.01
    OPEN
    SLIP
    warning('Slip is less than 0.01 in absolute value\n');
end

E =     xyobs(1,:) - pg(6);               % relative position of obs point wrt to fault centroid
N =     xyobs(2,:) - pg(7);               % relative position of obs point wrt to fault centroid

Emax = max(abs(E));
Nmax = max(abs(N));
if Emax > 100.e3 || Nmax > 100.e3 || hypot(Emax,Nmax) > 100.e3
    warning('Observation point is more than 100 km from fault centroid:');
    Emax
    Nmax
end



DEPTH = pg(3);                            % relative position of obs point wrt to fault centroid

NU    = nu;                               % Poisson's ratio

 
%
%	returns the following variables (same matrix size as E and N):
%	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
%	   uZE,uZN         : tilts (in rad * FACTOR)
%	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
 


%[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
%[uE,uN,uZ] = okada85(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU);
[uE,uN,uZ] = okada85disp(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU);

% figure
% hist(colvec(sqrt(uE.*uE + uN.*uN + uZ.*uZ)));
% title('displacement magnitude in meters');
% xlabel('U [m]');
% ylabel('number of occurrences');
% 
% disp 'uE'; size(uE)

uENZ(1,:) = uE;
uENZ(2,:) = uN;
uENZ(3,:) = uZ;

return



