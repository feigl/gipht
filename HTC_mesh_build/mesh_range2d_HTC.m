%% build a mesh from range change field
% 20170509 Kurt Feigl
close all
clear all
nf=0;

diary('run_HTC.log')

% get pair-specific variables from text file
% mast_yr mast_mn mast_day slav_yr slav_mn slav_day mpercy u_e u_n u_u
[VALS] = readtable('pair_info.txt', 'delimiter', ' ');

dt = dyear(table2array(VALS(1,4)),table2array(VALS(1,5)),table2array(VALS(1,6)))-dyear(table2array(VALS(1,1)),table2array(VALS(1,2)),table2array(VALS(1,3)))
depthsz =  50 % temperature

for dk = 1:numel(depthsz)
if fexist(sprintf('mesh_ranged2d%3dm.mat', -1*depthsz(dk))) == 0
    %% read range change
%     fn0 = 'drhomaskd_utm.grd';
    fn0 = 'x_001_OBSV_NUIS4.grd';
    INFO = grdinfo3(fn0)
    [datx,daty,grdr] = grdread3(fn0);
    grdr = grdr./1000; % put in mm
    ngrd = numel(grdr)
    
    %% assume grdr are in meters
    mpercy = table2array(VALS(1,7)); % meters per cycle
    
    
    %% make plaid grid with same spacing as interferogram
    [grdx,grdy]=meshgrid(datx,daty);
    vecx = colvec(grdx);
    vecy = colvec(grdy);
    vecr = colvec(grdr);
    
    %% select cells and prune
    ifast=find(vecr > 0.004);
    ifast=intersect(ifast,find(vecx >  327.0e3));
    ifast=intersect(ifast,find(vecy > 4405.5e3));
    fastx=vecx(ifast);
    fasty=vecy(ifast);
    fastr=vecr(ifast);
    
    %% unit vector pointing from ground to satellite
    % taken from /grdr/bradys/TSX/raw/20150714_orb.txt
    % pointing vector
    unitv_east   =  table2array(VALS(1,8))
    unitv_north  =  table2array(VALS(1,9))
    unitv_up     =  table2array(VALS(1,10))

    
    
    %% make plaid mesh with specified spacing in meters
    dx = 100;
    dy = 100;
    dz = 100;
    
    %% initial volume
    V0 = unit(dx*dy*dz,'m^3')
    
    gx = [min(vecx):dx:max(vecx)]';
    gy = [min(vecy):dy:max(vecy)]';
    % top and bottom elevations with respect to land surface
    topz = -1*depthsz(dk) % stay just below land surfact to avoid numerical instabilities
    botz = topz - 1 * dz;
    gz = (botz + topz)/2.;
    nz = numel(gz)
    ny = numel(gy)
    nx = numel(gx)
    [GX,GY] = meshgrid(gx,gy);
    % depth (currently constant 100)
    GZ = repmat(gz,size(GX));
    
    
%     %% make convex hull in 2 dimensions
%     %  make Delaunay triangulatiion
%     DelTriang = delaunayTriangulation(fastx,fasty);
%     % get extorior face of the region
%     khull = convexHull(DelTriang);
%     hullx = fastx(khull);
%     hully = fasty(khull);
    
    %% initialize constants
    V0 = unit(dx * dy * dz,'m^3');      % initial volume in cubic meters
    nu = 0.25               % Poisson's ratio
    vfactor = -1 %-1.e-6;          % scale factor for partial derivatives
    
    
    %% select locations to parameterize sinks
    iin = 1:numel(GZ);
    ncells = numel(iin)
    
    %% set centroids
    centx = colvec(GX(iin));
    centy = colvec(GY(iin));
    centz = colvec(GZ(iin));
    
%     %% build Gradient  operator
%     [DEL] = gradient_triangulate2d(centx,centy);

    %% build design matrix
    G = nan(ngrd,ncells);
    T=table(zeros(ncells,1),zeros(ncells,1),zeros(ncells,1),zeros(ncells,1),zeros(ncells,1),zeros(ncells,1),zeros(ncells,1) ...
        ,'VariableNames',{'x','y','z','dx','dy','dz','dV'});
    F=table2struct(T)
    %% calculate partial derivatives
    for i=1:ncells
        fprintf(1,'Cell number %d of %d\n',i,ncells);
        % axis are [easting, northing, elevation]
        % F.x,F.y,F.z cooardinates of fault centroids in meters (easting, northing, elevation) positive
        F(i).x = centx(i);
        F(i).y = centy(i);
        F(i).z = centz(i);
        % F(i).dx,F(i).dy,F(i).dz % dimensions of deforming rectangular prism
        F(i).dx = dx;
        F(i).dy = dy;
        F(i).dz = dz;
        % increasing volume is positive
        F(i).dV = vfactor * V0.value;
        % displacement 
        uENZ = okada85_wrapper3(rowvec(vecx),rowvec(vecy),F(i),nu); % use Kurt's function
        % okada displacement dotted with look vector of satellite
        dr = -1*(unitv_east * uENZ(1,:) + unitv_north * uENZ(2,:) + unitv_up * uENZ(3,:));
        G(:,i) = colvec(dr);
    end
    
    
    %% save entire structure
   % save('mesh_range2d.mat');
    save(sprintf('mesh_ranged2d%3dm.mat', topz));
    
end

%% load
topz=-1*depthsz(dk)
load(sprintf('mesh_ranged2d%3dm.mat', topz));

% fn0 = 'drhomaskd_utm.grd';
fn0 = 'x_001_OBSV_NUIS4.grd';
% fn0 = 'drhomaskd_utm.grd';
INFO = grdinfo3(fn0)
[datx,daty,grdr] = grdread3(fn0);
grdr=double(grdr);
grdr = grdr./1000; % put in mm
ngrd1 = numel(grdr);
if ngrd1 ~= ngrd
    nrgd1
    nrgrd
    error('miscount');
else
    ngrd=ngrd1
end

%% assume grdr are in meters
mpercy = table2array(VALS(1,7)); % meters per cycle

%% make plaid grid with same spacing as interferogram
[grdx,grdy]=meshgrid(datx,daty);
vecx = colvec(grdx);
vecy = colvec(grdy);
vecr = colvec(grdr);


%% prune missing data
ikeep = find(isfinite(vecr) == 1);

%% find a patch for reference
% refx =  327.e3; % UTM easting in meters
% refy = 4409.e3; % UTM northing in meters
% irefcount = 3000;
refx =  328.e3; % UTM easting in meters
refy = 4407.5e3; % UTM northing in meters
irefcount = 1000;
iref = find(hypot(vecx-refx,vecy-refy) < irefcount); %find(hypot(vecx-refx,vecy-refy) < 500); %find(hypot(vecx-refx,vecy-refy) < 3000);
iref = intersect(iref,ikeep);
refmed = nanmean(vecr(iref))

%% prune
drho_obs = vecr(ikeep);
datx_vec = vecx(ikeep);     
daty_vec = vecy(ikeep); 
G = G(ikeep,:);

A = G;
clear G;

drho_obs = drho_obs./dt;
b = drho_obs;
% 
% load('vecel.mat');
% vecel = vecel(ikeep);
% 
% mest = lscov([A, vecel], b);
% 
% b = b - 1.5*vecel*mest(end);

%% get deforming region
defregionx = [326.7e3, 327.8e3, 328.75e3, 329.25e3, 327.5e3, 326.7e3]';
defregiony = [4406.2e3, 4408.75e3, 4409.25e3, 4408.8e3, 4405.75e3, 4406.1e3]';
[iref, irefon] = inpolygon(centx, centy, defregionx, defregiony);
khull = find(iref > 0);
notkhull = find(iref == 0);
khullon = find(irefon > 0);
khull = [khull; khullon];

%% Geostatistical approach
% Reference Mike Cardiff's notes on geostatistical inverison
% X matrix for beta (assume no trend)
% Xb = ones(ncells, 1);
Xb = zeros(ncells, 2);
Xb(khull, 1) = 1;
Xb(notkhull, 2) = 1;

% set ranges for x and y
lx = 200 %500 %m, pressure prior
ly = 1000 %500 %m, pressure prior

% Perform Geostatistical Inversion
% unconstrained solution
[nrowsA,mcolsA] = size(A);
[nrowsb,mcolsb] = size(b);

[mdum, zcol] = size(A*Xb);

% save to mat files
save('A.mat', 'A');
save('b.mat', 'b');
save('Xb.mat', 'Xb');

% loop between temperature and pressure priors
for bi = 1
% define model covariance
if bi == 1
% temperature prior
disp('temperature prior')
mu = -4.2000e-05
sig2=(2.9086e-05)^2

else
%% pressure prior
disp('pressure prior')
mu=-9.1653e-06
sig2=(4.8109e-06)^2
end
Q = zeros(size(ncells));

% rotate into porotomo coordinate system
[centxp, centyp] = utm2xy_porotomo(centx, centy);
% 
% for i = 1:(ncells)
%     for j =1:i
%         Q(i,j) = sig2*exp(-sqrt((abs(centxp(i)-centxp(j))/(lx/3))^2 + (abs(centyp(i)-centyp(j))/(ly/3))^2));
%         Q(j,i) = Q(i,j);
%     end
% end

sig_nodef = (1e-20);

for i = 1:(ncells)
    for j =1:i
      if (numel(find(khull == (i))) > 0 && numel(find(khull == (j))) == 0) || numel(find(khull == (j))) == 0
        Q(i,j) = sig_nodef*exp(-sqrt((abs(centxp(i)-centxp(j))/(lx/3))^2 + (abs(centyp(i)-centyp(j))/(ly/3))^2));
        Q(j,i) = Q(i,j);
      else
        Q(i,j) = sig2*exp(-sqrt((abs(centxp(i)-centxp(j))/(lx/3))^2 + (abs(centyp(i)-centyp(j))/(ly/3))^2));
        Q(j,i) = Q(i,j);
%         int_ind(end+1) = i;
      end
    end
end

save('Q.mat', 'Q');

sizeb = numel(b);

data_unc  = .005; % uncertainty based on GPS

% clear workspace for room
clearvars -except A V0 vfactor zcol Q Xb nf depthsz centx centy dx dy datx daty nrowsA ncolsA datx_vec daty_vec grdx grdy dz dt sizeb irefcount data_unc sig2 mu bi ncells lx ly irefcount zcol topz

%% build A'*R where R is data spatial covariance matrix (too large to store)
% true representation
include_AR = 1
if include_AR == 1
%% build A'*R where R is data spatial covariance matrix (too large to store)
if fexist(sprintf('AR%3d.mat', topz)) == 0
c = 0.01*(0.005/dt)^2; %3e-6; %5e-8;
a = 230; %m
[arrows, arcols] = size(A');
arcols = sizeb;
AR = zeros(size(A'));
disp('Start defining AR')
start_time = clock
% for arj = 1:arcols %ari = 1:arrows
%     datacov_col = data_unc^(2)*(c - (c.*(1-exp(-3.*(abs(sqrt((datx_vec(arj) - datx_vec(:)).^2+(daty_vec(arj) - daty_vec(:)).^2))./a)))));
%     % define column of data covariance using spatial semivariogram
%     % model first
%     AR(:, arj) = A'*datacov_col;
% end
rind = [(0:3000:arcols), arcols];
for i = 1:numel(rind)-1
   AR(:, rind(i)+1:rind(i+1)) = A'*[c*exp(-3.*(abs(sqrt((repmat(datx_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(datx_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2+(repmat(daty_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(daty_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2))./a))];
end
ARfac = 1;
AR = AR*ARfac;
end_time = etime(clock, start_time)
clear datacov_col
save(sprintf('AR%3d.mat', topz), 'AR')
else
load(sprintf('AR%3d.mat', topz));
end
end

disp('Start defining Ggeo and Dgeo')

if include_AR == 0
% for no data covariance
% multiply by A and Xb to decrease size
 Ggeo = [(A'*A*Q)*A'+Xb*(A*Xb)', A'*(A*Xb)+Xb*(zeros(size(zcol)))];
else
% for when we have data covariance
% multiply by A and Xb to decrease size
C = zeros(zcol);
C(:,2) = 1; % constrain zero mean to non def region
Ggeo = [(A'*A*Q)*A'+AR+Xb*(A*Xb)', A'*(A*Xb)+Xb*C];
%clear AR
end
clear Q
disp('Ggeo defined')

load('b.mat')
Dgeo = A'*b+Xb*(zeros(zcol, 1));
clear b
disp('Dgeo defined')

clearvars -except A AR ARfac Ggeo Dgeo include_AR V0 vfactor nf depthsz centx centy dx dy datx daty nrowsA ncolsA datx_vec daty_vec grdx grdy dz dt sig2 mu bi ncells lx ly irefcount zcol topz


 % from regular ls with full dVdt3=unit(-23107,'m^3/year'), from 4Okada -3.1624e+04
 % from regular ls with < 500, DV =  -37107.3 m^3
% from Ali et al. (2016) result for In20130513_20140511: 2013.3616 2014.3562 -19244.912 3863.12 TSX T53_32785_38296 0.898067 0.906401 0.9946
%pest = pinv(Ggeo'*Ggeo)*Ggeo'*Dgeo; %doulbe version: DV = -376305 m^3 %pinv(Ggeo'*Ggeo,  1e-18)*Ggeo'*Dgeo; % DV =  -37236.9 m^3
load('Q.mat')
disp('Start inversion')
pest = pinv((Ggeo))*(Dgeo);

clear Dgeo

disp('inversion completed')
load('A.mat')
load('Xb.mat')
load('Q.mat')

disp('finding posterior model covariance')

% store only diagonal elements
psig = diag(Q - [A*Q; Xb']'*pinv(Ggeo)*[A', Xb]*[A*Q; Xb']);

disp('defining beta and xi')
clear Ggeo
load('b.mat')
% beta = pest(end)
% xi = pest(1:end-1);
beta = pest(end-1:end)
xi = pest(1:end-2);

% mode 
disp('Mode of estimated parameter values')
mest = Xb*beta + Q*A'*xi;

% compare to true model parameters
if include_AR == 0
mse = (b-A*mest)'*(b-A*mest)
else
mse0 = (b-A*mest)'*(b-A*mest)
res = (b-A*mest);
load('AR-50.mat')
[Arows, Acols] = size(A);
degfree = Arows - Acols;
mse_pre = ((res)'*(res))/((.005/dt)^2*ARfac)/degfree
mse = (b-A*mest)'*pinv(AR)*A'*(b-A*mest)*ARfac
clear AR
end
irun = 1;

% find likelihood and prior for bayes factor
if bi == 1
llikelihood_dT = -1/2*mse
lprior_dT = -1/2*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)
lpost_dT = lprior_dT*llikelihood_dT
end

if bi == 2
llikelihood_dP = -1/2*mse
lprior_dP = -1/2*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)
lpost_dP = lprior_dP*llikelihood_dP
end



%% Material properties
nu = 0.2;  % dimensionless Poisson's ratio
dVdT = unit(30.e-6,'1/K') ; % volumetric coefficient of thermal expansion in degrees per Kelvin
G_shear = unit(3.e9,'Pa') % in Pa

% calculate Bulk Modulus
K_bulk_mod = 2*G_shear*(1+nu)/(3*(1-nu)) % 3.3333e+10 Pa
% Nur and Byerlee [1971] give K = 36.E5 bar = 3.6E10 Pa for quartz, citing Birch [1966]
% Ingraham et al. give a median value of 8466.5 MPa = 8.4E9 Pa

% Biot-Willis coefficient K/H Wang (1.24), page 21
%Ingraham, M. D., S. J. Bauer, K. A. Issen, and T. A. Dewers (2017),
%Evolution of permeability and Biot coefficient at high mean stresses in
%high porosity sandstone, International Journal of Rock Mechanics and
%Mining Sciences, 96, 1-10. They give a median value of:
alpha_biot = 0.8505

% Poroelastic expansion coefficient 1/H
invH = alpha_biot / K_bulk_mod


%% Change in Volume
mest = double(mest);
disp('vfactor = ')
vfactor
dV = vfactor * mest * V0;
dV_std = vfactor * sqrt(abs(psig)) * V0;
disp('V0 = ')
V0

%% total
DV = unit(sum(dV.value),'m^3')
DV_defregion = sum(dV.value((dV.value) < -2*std(dV.value)))
DV_defregion_std = sum(dV_std.value((dV.value) < -2*std(dV.value)))

%% total from reservoir region TSX T53 In20160722_20170822
load('resv_ind.mat')
DV_resregion = sum(dV.value(resv_ind))
DV_resregion_std = sum(dV_std.value(resv_ind))

%% temperature change
dT = dV / V0 / dVdT;
dT_std = dV_std / V0 / dVdT;
if bi == 1
disp('estimated product of dT*eps')
dTe = mean(mest)
dTe_std = std(mest)
fprintf('bounds = (%5.5e, %5.5e) \n', mu-sqrt(sig2), mu+sqrt(sig2))
end

dTmean = convert(unit(nanmean(dT.value),'degC'),'degC')
dTmean_defregion = mean(dT.value(abs(dT.value) > 2*std(dT.value)))
DTmean_defregion_std = sum(dT_std.value(abs(dT.value) > 2*std(dT.value)))

%% pressure change
dP = dV / V0 / (invH );
if bi == 2
disp('estimated product of dP*1/H')
%dPH = mean(dV.value/ V0.value)
dPH = mean(mest)
%dPH_std = std(dV.value/V0.value)
dPH_std = std(mest)
fprintf('bounds = (%5.5e, %5.5e) \n', mu-sqrt(sig2), mu+sqrt(sig2))
end
dPmean = unit(nanmean(dP.value),'Pa')


%% stress change
dS = K_bulk_mod.value * dV/V0;
dSmean = unit(nanmean(dS.value),'Pa')

%% strain rate
dStrainmean_defregion = mean(dV.value((dV.value) < -2*std(dV.value)))/V0/3.154e7
dStrainmean_defregion_std = mse0*dStrainmean_defregion

%% calculate Energy rate 
Cs = unit(948, 'J/degC/kg'); % specific heat of rock, Rutqvist
%Cs = unit(800, 'J/degC/kg'); % specific heat of rock, Ali, et al (2016)
mparams = numel(mest);
drho_lb = unit(1900, 'kg/m^3') ; % kg/m^3 lower bound for rock density from Witter
drho_ub = unit(2800, 'kg/m^3'); % kg/m^3 upper bound for rock density from Witter
% dT_std = std(dT.value)
% dT2sig = find(dT.value > 2*dT_std);
% dT2sigtotal = sum(dT.value(dT2sig))
% dTtotal = sum(dT.value) % total temperature change rate

%dEdt_lb = Cs*mparams*V0*drho_lb*dTtotal
%dEdt_ub = Cs*mparams*V0*drho_ub*dTtotal
sumdV = unit(sum(dV.value((dV.value) < -2*std(dV.value))), 'm^3/year')
meandV = unit(mean(dV.value((dV.value) < -2*std(dV.value))), 'm^3/year')
stddV =  unit(sum(dV_std.value((dV.value) < -2*std(dV.value))), 'm^3/year')

% 1 watt = 3.16887646e8 J/yr
%dE_lb = unit(Cs.value*1/dVdT.value*drho_lb.value*sum(dV.value)*3.2e-8, 'watts')
%dE_ub = unit(Cs.value*1/dVdT.value*drho_ub.value*sum(dV.value)*3.2e-8, 'watts')
dE_lb = Cs*1/dVdT*drho_lb*sumdV;
dE_ub = Cs*1/dVdT*drho_ub*sumdV;
dE_mean = Cs*1/dVdT*drho_ub*meandV;
dE_std = Cs*1/dVdT*drho_ub*stddV;
dE_lb = convert(dE_lb, 'watts')
dE_ub = convert(dE_ub, 'watts')
dE_mean = convert(dE_mean, 'watts')
dE_std = convert(dE_std, 'watts')


%%for shorter run without plots end here
if bi == 1
    save('mest_dT.mat', 'mest', 'psig')
%     dlmwrite('dV_dT.txt', [dV.value, dV_std.value]) 
    dlmwrite('mest_dT.txt', [mest, psig], 'delimiter', ' ') 
else
    save('mest_dP.mat', 'mest')
    dlmwrite('dV_dP.txt', [dV.value, dV_std.value]) 
end
% end loop over priors
end
end
diary off









