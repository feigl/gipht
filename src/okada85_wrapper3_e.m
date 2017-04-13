function uENZ = okada85_wrapper3_e(F, xyobs)
%function uENZ = okada85_wrapper3(xobs,yobs,F,nu)
% wrapper around Okada85 routines for three orthogonal dikes

%fprintf(1,'entering %s at %s\n',mfilename,datestr(now));
%% inputs
% axis are [easting, northing, elevation]
% F.xcentroid,F.ycentroid,F.zcentroid cooardinates of fault centroids in meters 
% F.dx,F.dy,F.dz % dimensions of deforming rectangular prism
% F.vstrain      % dimensionless strain: increasing volume is positive
xobs = xyobs(:,1);
yobs = xyobs(:,2);
%% count number of observation points
nobs= numel(xobs);

%% count number of fault centroids
ncentroids = numel(F.xcentroid);

%% initialize
uENZ = zeros(3,nobs);

%% common fault parameters
rake    =    0; % rake of slip vector in degrees
slip    =    0; % amount of in-plane slip in meters
%nu           ; % Poisson's ratio


% whos
% F

%% change in volume [m^3]
% F.volchange = F.volstrain .* (F.dx .* F.dy .* F.dz);
F.volchange = F.volstrain;

%% loop over fault centroids
onev = ones(size(F.xcentroid));
%for j=1:ncentroids
for j=1
    %% coordinates of observation points in meters
    % make observation point on ground above
    Eobspts_wrt_fault =     xobs  - F.xcentroid(j);               % relative position of obs point wrt to fault centroid
    Nobspts_wrt_fault =     yobs - F.ycentroid(j);               % relative position of obs point wrt to fault centroid
    
    %% depth of fault centroid in meters
    depth     =    -1*F.zcentroid(j) ;                   % vertical coordinate of fault centroid   
    
    %% calculate opening from volume strain
    %volume = F.dx(j) .* F.dy(j) .* F.dz(j);          % cubic meters
    
    %% loop over three orthogonal dikes
    for i=1:3
        switch i
            case 1 % horizontal fault
                dipd    =    F.dip+0;  % dip in degrees
                striked =    F.strike+0;  % strike in degrees
                length  = F.dx(j);  % length of fault in meters
                width   = F.dy(j);  % width of fault in meters
            case 2 % vertical fault in X-Z plane
                 dipd    =   F.dip+90;  % dip in degrees
                striked =   F.strike+90;  % strike in degrees clockwise from north
                length  = F.dx(j); % length of fault in meters
                width   = F.dz(j); % width of fault in meters
            case 3 % vertical fault in X-Y plane
                dipd    =   F.dip+90;  % dip in degrees
                striked =    F.strike+0;  % strike in degrees clockwise from north
                length  = F.dy(j); % length of fault in meters
                width   = F.dz(j); % width of fault in meters
            otherwise
                error(sprintf('unknown value of i = %d\n',i));
        end
        if abs(depth) < width
            error('depth < width')
        end
        %% tensile opening on one dike
        open   = F.volchange(j) / 3.0 / length / width; % meters
        
%         %% make vectors
%         depth   = depth   * onev;
%         dipd    = dipd    * onev;
%         striked = striked * onev;
%         length  = length  * onev;
%         width   = width   * onev;
%         rake    = rake    * onev;
%         slip    = slip    * onev;
%         open    = open    * onev;
        
        
        %% run calculate displacements from one set of dikes
        [uE,uN,uZ] = okada85disp(Eobspts_wrt_fault,Nobspts_wrt_fault,depth,striked,dipd,length,width,rake,slip,open,F.nu);
        
        %% sum the displacements over three dikes
        uENZ(1,:) = uENZ(1,:) + uE';
        uENZ(2,:) = uENZ(2,:) + uN';
        uENZ(3,:) = uENZ(3,:) + uZ';
   end
    
end

%whos
return



