function uENZ = okada85_wrapperlwh_brady(xobs,yobs,F,nu)
% computes a brick/book shape with different l, w, and h dimensions with
% each plane modeled by an Okada dislocation source
% Elena C Reinisch 20180827

%fprintf(1,'entering %s at %s\n',mfilename,datestr(now));
%% inputs
% axis are [easting, northing, elevation]
% F(j).x,F(j).y,F(j).z cooardinates of fault centroids in meters (easting, northing, elevation) positive 
% F(j).dx,F(j).dy,F(j).dz % dimensions of deforming rectangular prism
% F(j).volume_increase      % dimensionless strain: increasing volume is positive

%% count number of observation points
nobs= numel(xobs);

%% count number of fault centroids
ncentroids = numel(F);

%% initialize
uENZ = zeros(3,nobs);

%% common fault parameters
rake    =    0; % rake of slip vector in degrees
slip    =    0; % amount of in-plane slip in meters
%nu           ; % Poisson's ratio


% whos
% F


%% loop over fault centroids
%onev = ones(size(F(j).x));
% for j=1
for j=1:ncentroids
    %% coordinates of observation points in meters
    Eobspts_wrt_fault =     xobs - F(j).x;               % relative position of obs point wrt to fault centroid
    Nobspts_wrt_fault =     yobs - F(j).y;               % relative position of obs point wrt to fault centroid
    
    %% depth (below surface) of fault centroid in meters
    %depth     =    F(j).z ;                  
    % 20170207 corrected so that F(j).z is elevation above ground surface
    % 20170509 confirmed Kurt
    depth     =    -1*F(j).z ;                  
    
    %% calculate opening from volume strain
    %volume = F(j).dx .* F(j).dy .* F(j).dz;          % cubic meters
    %% calculate dimensions of volume change in terms of length assuming ratio of sides stays the same
    % dV = V2 - V0
    % dV = L0*W0*H0 - L2*W2*H2
    V0 = F(j).dx*F(j).dy*F(j).dz;
    L2 = ((F(j).dV + V0)/((F(j).dy/F(j).dx)*(F(j).dz/F(j).dx)))^(1/3);
    W2 = (F(j).dy/F(j).dx)*L2;
    H2 = (F(j).dz/F(j).dx)*L2;
    
    
    %% loop over three orthogonal dikes
    for i=1:3
        switch i
            case 1 % horizontal fault
                dipd    =    0;  % dip in degrees
                striked =    36.4219;  % strike in degrees
                length  = F(j).dx;  % length of fault in meters
                width   = F(j).dy;  % width of fault in meters
                open   = H2 - F(j).dz;
            case 2 % vertical fault in X-Z plane
                dipd    =   90;  % dip in degrees
                striked =   36.4219+90;  % strike in degrees clockwise from north
                length  = F(j).dx; % length of fault in meters
                width   = F(j).dz; % width of fault in meters
                open   = W2 - F(j).dy;
            case 3 % vertical fault in X-Y plane
                dipd    =   90;  % dip in degrees
                striked =    36.4219;  % strike in degrees clockwise from north
                length  = F(j).dy; % length of fault in meters
                width   = F(j).dz; % width of fault in meters
                open   = L2 - F(j).dx;
            otherwise
                error(sprintf('unknown value of i = %d\n',i));
        end
        
        %% tensile opening on one dike
%         open   = F(j).dV / 3.0 / length / width; % meters
        
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
        [uE,uN,uZ] = okada85disp(Eobspts_wrt_fault,Nobspts_wrt_fault,depth,striked,dipd,length,width,rake,slip,open,nu);
        
        %% sum the displacements over three dikes
        uENZ(1,:) = uENZ(1,:) + uE;
        uENZ(2,:) = uENZ(2,:) + uN;
        uENZ(3,:) = uENZ(3,:) + uZ;
    end
    
end

%whos
return



