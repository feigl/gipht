function DISPLACEMENT = run_okada_for_volume_change(SSS,t0,t1,make_plots)

%% run the Okada Model
% 20170404 Kurt Feigl

%% initialize
tcomp = tic;
nf = 0;


% 
%% Material properties
%SSS.nu = 0.25;  % dimensionless Poisson's ratio
%SSS.dVdT = 30.e-6 ; % volumetric coefficient of thermal expansion in degrees per Kelvin
%SSS.G_shear = 3.e10 % in Pa

%% calculate Bulk Modulus 
K_bulk_mod = 2*SSS.G_shear*(1+SSS.nu)/(3*(1-SSS.nu; % 3.3333e+10 Pa 
% Nur and Byerlee [1971] give K = 36.E5 bar = 3.6E10 Pa for quartz, citing Birch [1966] 

%% check grid dimensions
[nvoxelsG,dummy] = size(SSS.GRID.GridBlock);
if numel(SSS.GRID.CornersX) ~= 8*nvoxelsG
    error('Miscount in X');
end
if numel(SSS.GRID.CornersY) ~= 8*nvoxelsG
    error('Miscount in Y');
end
if numel(SSS.GRID.CornersZ) ~= 8*nvoxelsG
    error('Miscount in Z');
end


%% coordinates are model coordinates in meters
% xobs = linspace(-5e3,5e3,100);
% yobs = linspace(-5e3,5e3,100);
xaxe = linspace(SSS.GRID.xmin,SSS.GRID.xmax,200);
yaxe = linspace(SSS.GRID.ymin,SSS.GRID.ymax,240);
[xobs2,yobs2] = meshgrid(xaxe,yaxe);
xobs = rowvec(xobs2);
yobs = rowvec(yobs2);

%% top
zmin = min(SSS.GRID.ModelZ);

%% make a plot with grid points
if make_plots == 1
    figure;
    plot(xobs,yobs,'k.');
    xlabel('Xmodel [m]');
    ylabel('Ymodel [m]');
    title('xobs and yobs');
end

%% calculate rate of change of Temperature in degC per year
dyears = years(diff(SSS.PTFIELD.t'));
dTdt = diff(SSS.PTFIELD.T,1,2)./repmat(dyears,size(SSS.PTFIELD.T,1),1);
dTdt_mean = mean(dTdt,2);

%% calculate rate of change of Pressure in Pascal per year
dyears = years(diff(SSS.PTFIELD.t'));
dPdt = diff(SSS.PTFIELD.P,1,2)./repmat(dyears,size(SSS.PTFIELD.P,1),1);
dPdt_mean = mean(dPdt,2);




%% select voxels (gridblocks)
ii = 1:nvoxelsG; % all voxels
%ii = ifast_voxels;
%ii = 1:100; %shorten run for debugging
ii = intersect(ii,find(abs(dTdt_mean) > SSS.Tthresh)); % rapidly changing pixels
%ii = intersect(ii,find(dTdt_mean < -0.5)); % contracting only
ii = intersect(ii,find(abs(dPdt_mean) > SSS.Pthresh)); % rapidly changing pixels

nvoxels = numel(ii);

if make_plots == 1
    %% histogram of temperature change
    nf=nf+1;h(nf)=figure;hold on;
    histogram(dTdt_mean(ii),30);
    xlabel('dTdt [degC/year]');
    ylabel('number of voxels');
    title('Rate of change of Temperature Change');
    printpdf(sprintf('%s_%02d.pdf',mfilename,nf));
    
    %% histogram of pressure change
    nf=nf+1;h(nf)=figure;hold on;
    histogram(dPdt_mean(ii),30);
    xlabel('dPdt [Pa/year]');
    ylabel('number of voxels');
    title('Rate of change of Pressure Change');
    printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

end


%% calculate displacement vector in meters for each voxel
u_xyz = zeros(3,numel(xobs));
tfirst = datetime(0,1,1);
tfirst.TimeZone='UTC';
tfirst.Format='yyyyMMdd';
for i=1:nvoxels
    % coordinates of centroids of voxels
    F.x = SSS.GRID.ModelX(ii(i));
    F.y = SSS.GRID.ModelY(ii(i));
    %F.z = SSS.GRID.ModelZ(ii(i)) - zmin; % depth below surface, positive downward
    F.z = zmin - SSS.GRID.ModelZ(ii(i)); % elevation above surface, positive upward
    % dimensions of voxels
    F.dx = SSS.GRID.ModelDX_m_(ii(i));
    F.dy = SSS.GRID.ModelDY_m_(ii(i));
    F.dz = SSS.GRID.ModelDZ_m_(ii(i));
    
    %% original volume
    V0 = F.dx .* F.dy .* F.dz;
    
    %% change in volume [m^3] from thermal expansion
    if SSS.dVdT > 0.       
        %% change in temperature in degrees Celsius
%         T0 = interp1(years(SSS.PTFIELD.t-datetime(0,1,1)),SSS.PTFIELD.T(ii(i),:)',years(SSS.t0-datetime(0,1,1)));
%         T1 = interp1(years(SSS.PTFIELD.t-datetime(0,1,1)),SSS.PTFIELD.T(ii(i),:)',years(SSS.t1-datetime(0,1,1)));
        T0 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.T(ii(i),:)',years(t0-tfirst));
        T1 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.T(ii(i),:)',years(t1-tfirst));
        dT = T1 - T0;
        
        %% volumetric strain [dimensionless]
        volume_strainT = dT * SSS.dVdT;  
    else
        dT = 0;
        volume_strainT = zeros(size(V0));  
    end
    
    %% change in volume [m^3] from thermal expansion
    %F.volume_increase = volume_strain * V0;
    
    %     %% try Okada solution
    %     u_xyz = u_xyz + okada85_wrapper3(xobs,yobs,F,SSS.nu);
    
    %% change in volume [m^3] from pressure change
    if SSS.G_shear > 0. && isfinite(K_bulk_mod) == 1
%         P0 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.P(ii(i),:)',years(SSS.t0-tfirst));
%         P1 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.P(ii(i),:)',years(SSS.t1-tfirst));
        P0 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.P(ii(i),:)',years(t0-tfirst));
        P1 = interp1(years(SSS.PTFIELD.t-tfirst),SSS.PTFIELD.P(ii(i),:)',years(t1-tfirst));
        dP = P1 - P0;
        
        %% volumetric strain [dimensionless]
        % Taking increasing volume as positive strain, increasing confining
        % pressure on outside of a unit cube should DECREASE volume. This
        % gives subsidence at Brady
        % volume_strainP = -1 * dP / K_bulk_mod;
        % 20170403 Taking increasing volume as positive strain, increasing
        % pressure of fluids in pores should increase volume. This gives
        % uplift at Brady
        volume_strainP = dP / K_bulk_mod;
    else
        dP = 0;
        volume_strainP = zeros(size(V0));      
    end
    
    %% change in volume [m^3] from pressure and volume
    F.volume_increase = (volume_strainT + volume_strainP) * V0;
    
    %% run the Okada solution - depends on Poisson ratio but not shear modulus
    u_xyz = u_xyz + okada85_wrapper3(xobs,yobs,F,SSS.nu);
     
%     fprintf(1,'i = %6d voxel ii(i) = %6d depth %6.1f m dT = %7.3f degC dP = %10.2E Pa dV = %10.2E m^3\n',i,ii(i),F.dz,dT,dP,F.volume_increase);
   
    %% let us know how we are doing
    if mod(i,100) == 0
        fprintf(1,'Completed computing voxel %6d of %6d after %s\n',i,nvoxels,seconds(toc(tcomp)));
    end
end

%% now rotate modeled coordinates and displacements into UTM easting, northing
fprintf(1,'Calculating rotated coordinates for point number:');
xyz=zeros(3,1);
for i=1:numel(xobs)
    if mod(i,1000) == 0
        fprintf(1,' %d',i);
    end
    xyz(1)=xobs(i)-SSS.ORIGIN.ModelX;
    xyz(2)=yobs(i)-SSS.ORIGIN.ModelY;
    xyz(3)=0.;  % assume flat earth
    %% rotate coordinates
    enu=SSS.GRID.R3R2*xyz;
    DISPLACEMENT.e(i) = enu(1) + SSS.ORIGIN.Easting;
    DISPLACEMENT.n(i) = enu(2) + SSS.ORIGIN.Northing;
    DISPLACEMENT.v(i) = enu(3) + SSS.ORIGIN.Elevation;
    %% rotate displacements
    xyz(1)=u_xyz(1,i);
    xyz(2)=u_xyz(2,i);
    xyz(3)=u_xyz(3,i);
    enu=SSS.GRID.R3R2*xyz;
    DISPLACEMENT.ue(i) = enu(1);
    DISPLACEMENT.un(i) = enu(2);
   %DISPLACEMENT.uv(i) = enu(3); % identified as wrong 20170121
    DISPLACEMENT.uv(i) =-1*enu(3); % corrected 20170121
end
fprintf(1,'\nDone.\n');


%% archive the results
[DISPLACEMENT.nrows, DISPLACEMENT.ncols] = size(xobs2);
DISPLACEMENT.e = colvec(DISPLACEMENT.e);
DISPLACEMENT.n = colvec(DISPLACEMENT.n);
DISPLACEMENT.v = colvec(DISPLACEMENT.v);
DISPLACEMENT.ue = colvec(DISPLACEMENT.ue);
DISPLACEMENT.un = colvec(DISPLACEMENT.un);
DISPLACEMENT.uv = colvec(DISPLACEMENT.uv);
DISPLACEMENT.x  = colvec(xobs);
DISPLACEMENT.y  = colvec(yobs);
DISPLACEMENT.ux = colvec(u_xyz(1,:));
DISPLACEMENT.uy = colvec(u_xyz(2,:));
DISPLACEMENT.uz = colvec(u_xyz(3,:));

%save(sprintf('%s_DISPLACEMENT.mat',mfilename),'-struct','DISPLACEMENT')

return


