function uENZ = okada85_wrapper3_e(xobs,yobs,F,nu)
%function uENZ = okada85_wrapper3_e(xobs,yobs,F,nu)
% wrapper around Okada85 routines for three orthogonal dikes with output of
% regular surface displacement (non-InSAR convention, deltaZ < 0 indicates subsidence, see mogi_e.m for
% analogous solution)
% inputs
% axis are [easting, northing, elevation]
% F.x,F.y,F.z cooardinates of fault centroids in meters, z is positive
% depth
% F.dx,F.dy,F.dz % dimensions of deforming rectangular prism
% F.vstrain      % dimensionless strain: increasing volume is positive
%
% edited Elena C Reinisch 20170505 change to input of positive depth and
% output of surface displacement with deltaZ < 0 --> subsidence convention

%% count number of observation points
nobs= numel(xobs);

%% count number of fault centroids
ncentroids = numel(F.x);

%% initialize
uENZ = zeros(3,nobs);

%% common fault parameters
rake    =    0; % rake of slip vector in degrees
slip    =    0; % amount of in-plane slip in meters
%nu           ; % Poisson's ratio


% whos
% F


%% loop over fault centroids
onev = ones(size(F.x));
%for j=1:ncentroids
for j=1
    %% coordinates of observation points in meters
    Eobspts_wrt_fault =     xobs - F.x(j);               % relative position of obs point wrt to fault centroid
    Nobspts_wrt_fault =     yobs - F.y(j);               % relative position of obs point wrt to fault centroid
    
    %% depth (below surface) of fault centroid in meters
    depth     =    F.z(j) ;  % depth is positive below surface     
    
    %% calculate opening from volume strain
    %volume = F.dx(j) .* F.dy(j) .* F.dz(j);          % cubic meters
    
    %% loop over three orthogonal dikes
    for i=1:3
        switch i
            case 1 % horizontal fault
                dipd    =    0;  % dip in degrees
                striked =    0;  % strike in degrees
                length  = F.dx(j);  % length of fault in meters
                width   = F.dy(j);  % width of fault in meters
            case 2 % vertical fault in X-Z plane
                dipd    =   90;  % dip in degrees
                striked =   90;  % strike in degrees clockwise from north
                length  = F.dx(j); % length of fault in meters
                width   = F.dz(j); % width of fault in meters
            case 3 % vertical fault in X-Y plane
                dipd    =   90;  % dip in degrees
                striked =    0;  % strike in degrees clockwise from north
                length  = F.dy(j); % length of fault in meters
                width   = F.dz(j); % width of fault in meters
            otherwise
                error(sprintf('unknown value of i = %d\n',i));
        end
        
        %% tensile opening on one dike
        open   = F.volume_increase(j) / 3.0 / length / width; % meters
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
%         disp('SIGN DEPTH')
%         sign(depth)
        [uE,uN,uZ] = okada85disp(Eobspts_wrt_fault,Nobspts_wrt_fault,depth,striked,dipd,length,width,rake,slip,open,nu);

        %% multiply Z axis by -1 to change from positive depth convention to elevation (Bonafede, M., and C. Ferrari (2009), p 5)
        uZ = -1*uZ; 
        
        %% sum the displacements over three dikes
        uENZ(1,:) = uENZ(1,:) + uE;
        uENZ(2,:) = uENZ(2,:) + uN;
        uENZ(3,:) = uENZ(3,:) + uZ;
    end
    
end

%whos
return



