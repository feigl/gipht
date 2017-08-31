function uENZ = strain_nucleus(xobs,yobs,F,nu)
%function uENZ = strain_nucleus(xobs,yobs,F,nu)
% calculate displacement due to expansion of a strain nuclei
%% inputs
% coordinate axes are (xobs, yobs)
% vertical component of displacement is positive upwards
%  
% F.x,F.y,F.z cooardinates of fault centroids in meters
% F.dx,F.dy,F.dz % dimensions of deforming rectangular prism
% F.Vstrain      % volumetric strain DV/V [dimensionless]
% nu        % Poisson's ratio
% Kurt Feigl 20161006

%% count number of observation points
nobs= numel(xobs);

%% count number of fault centroids
nnuclei = numel(F.x);

%% initialize
uENZ = zeros(3,nobs);

nsteps = 10;

%% loop over observation points
for i=1:nobs
    
    %% loop over fault nuclei
    UHS = zeros(3,1);
    for j=1:nnuclei
        %% loop over integration points
        %% coordinates of integration points in meters
        Ix = linspace(F.x(j)-F.dx(j)/2., F.x(j)+F.dx(j)/2., nsteps);
        Iy = linspace(F.y(j)-F.dy(j)/2., F.y(j)+F.dy(j)/2., nsteps);
        Iz = linspace(F.z(j)-F.dz(j)/2., F.z(j)+F.dz(j)/2., nsteps); % internally modelZ is positive downward, i.e. depth]
        %dvolume_strain = F.volume_strain/(nsteps^3);
        dvolume_strain = F.volume_strain(j);
        
        %% sum over incremental volume in a single nucleus
        dUHS = zeros(3,1);
        for ix = 1:nsteps
            for iy = 1:nsteps
                for iz = 1:nsteps
                    Xx = xobs(i) - Ix(ix);  % relative position of observation point wrt to integration point in strain nucleus
                    Xy = yobs(i) - Iy(iy);  % relative position of observation point wrt to integration point in strain nucleus
                    zp =      0. - Iz(iz);  % relative position of observation point wrt to integration point in strain nucleus
                    % dUHS = dUHS + wang539(dvolume_strain,nu,zp,Xx,Xy);
                    dUHS = dUHS + wang539(F.dV(j)/(nsteps^3),nu,zp,Xx,Xy);
                end
            end
        end
        %% sum the displacements over nucleii
        
        UHS = UHS + dUHS;
    end
    
    %% pack the displacement vector for the observation point into array
    uENZ(1,i) =    UHS(1);
    uENZ(2,i) =    UHS(2);
    uENZ(3,i) = -1*UHS(3); % return vertical component positive updwards
    
end

return



