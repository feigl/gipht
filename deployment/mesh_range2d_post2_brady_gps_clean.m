%% build a mesh from GPS range change field and estimate displacements from cubiodal model
% last updated 20191013 Elena C Reinisch
close all
clear all

nf=0;

%% Get data
load('../GPS/BRD1_BRDY_new_full2.mat');
s_hat = [-0.5559;   0.1011;    -0.8251];
brady_range = GPS.BRD1.range - GPS.BRDY.range;
brady_range_unc = ([GPS.BRD1_BRDY.sig_east, GPS.BRD1_BRDY.sig_north, GPS.BRD1_BRDY.sig_vert])*s_hat/1e3;

% % form pairs
% View GPS dates
% convert dates to datetime format for comparison with InSAR
[tu_y, tu_m, tu_d] = dyear2date(GPS.BRD1_BRDY.tB1BY);

% select only pairs in deployment period
deployment_ind = 1:numel(tu_y);
tu_dt = datetime(tu_y(deployment_ind), tu_m(deployment_ind), tu_d(deployment_ind));
tu = GPS.BRD1_BRDY.tB1BY(deployment_ind);
brady_range = brady_range(deployment_ind);
brady_range_unc = brady_range_unc(deployment_ind);

% daily pairs
tm = tu(1:end-1);
ts = tu(2:end);
[tm_y, tm_m, tm_d] = dyear2yyyymmdd(tm);
tm_dt = datetime(tm_y, tm_m, tm_d);
[ts_y, ts_m, ts_d] = dyear2yyyymmdd(ts);
ts_dt = datetime(ts_y, ts_m, ts_d);
brady_drho = (brady_range(2:end) - brady_range(1:end-1))./(ts - tm);
brady_drho_unc = sqrt(brady_range_unc(2:end).^2 + brady_range_unc(1:end-1).^2)./(ts - tm);

brady_drho1_unc = sqrt(brady_range_unc(1:end-1).^2)./(ts - tm);
brady_drho2_unc = sqrt(brady_range_unc(2:end).^2)./(ts - tm);

ndat = numel(brady_drho);

depthsz =  50 

for dk = 1:numel(depthsz)
   [x_gps, y_gps] = deg2utm(39.809, -119.003);
    
    
    %% unit vector pointing from ground to satellite
    % pointing vector
    unitv_east   =  0.5559;   
    unitv_north  =  -0.1011;
    unitv_up     =  0.8251;

    
    
    %% make plaid mesh with specified spacing in meters
    dx = 1000;
    dy = 3000;
    dz = 100;
    topz = -1*depthsz(dk) % stay just below land surfact to avoid numerical instabilities

    
    %% initialize constants
    V0 = unit(dx * dy * dz,'m^3');      % initial volume in cubic meters
    nu = 2.5E-01              % Poisson's ratio
    vfactor = -1 %-1.e-6;          % scale factor for partial derivatives
    
   % set centroids
    centx = 3.279124094369082e+05;
    centy =  4.407489059350844e+06;
    centz = 100; 
    

    %% build design matrix
    ncells = 1;
    %G = nan(ndat,ncells);
    %G = nan(numel(grdx),ncells);
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
        uENZ = okada85_wrapperlwh_brady(x_gps, y_gps, F(i),nu); % use Kurt's function
        
        % okada displacement dotted with look vector of satellite
        dr = -1*(unitv_east * uENZ(1,:) + unitv_north * uENZ(2,:) + unitv_up * uENZ(3,:));
    end
    dV_est = -1*(brady_drho/dr)*V0.value;
    dV_unc = sqrt(abs((V0.value/dr).^2*brady_drho_unc.^2));
end        
    % save pairs
   dlmwrite('brady_vol_est_gps_ROI_new_full2.txt', [tm, ts, dV_est, dV_unc], 'delimiter', ' ', 'precision', '%.8e');
