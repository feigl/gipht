function [solution_epochs,xpts,ypts,zpts,dupts,dvpts,dwpts] = get_comsol_displacements(PST,DST,TST)
%function [solution_epochs,xpts,ypts,zpts,dupts,dvpts,dwpts] = get_comsol_displacements(PST,DST,TST)
% Kurt 2014-08-08
xpts=DST.x;
ypts=DST.y;
zpts=DST.z;

verbose = 0;
%verbose = 1;

%[pnames32, pnames, values, dims, descrs] = get_comsol_parameters(PST.datafilename)

% Comsol evaluation solution_epochs in seconds
%[model, solution_epochs] = execute_comsol(PST.datafilename)

model = mphload(PST.datafilename);

info1 = mphsolinfo(model);
info2 = mphsolutioninfo(model);

if verbose == 1
    info1
    info2
end

% Comsol evaluation solution_epochs in seconds
solution_epochs = info1.solvals;

% % set Poisson's ratio
% %PST.NUP=0.2;
% if isfinite(PST.NUP)==1
%     model.param.set('NUP', sprintf('%e',PST.NUP),'Poissons Ratio');
% end
%
% % set Young's Modulus
% if isfinite(PST.EYM)==1
%     model.param.set('EYM', sprintf('%e[Pa]',PST.EYM), 'Youngs Modulus');
% end
%
% model.param.set('LENG', sprintf('%e[m]',PST.LENG), 'ellipsoid a axis');
% model.param.set('WIDTH', sprintf('%e[m]',PST.WIDTH), 'ellipsoid b axis');
% %model.param.set('PRES', sprintf('%e[Pa]',PST.PRES), 'pressure');


%obtain nodal coordinates first:
nodestruct = mphxmeshinfo(model,'soltag','sol1');
coords = nodestruct.nodes.coords;
%max(max(coords));
[ndim,ncols] = size(coords);
x_node = coords(1,:); % x coordinate  in meters, w.r.t. xcen
y_node = coords(2,:); % y coordinate  in meters, w.r.t. ycen
z_node = coords(3,:); % z component in meters, positive upwards, w.r.t. 0
% size(x_node);
% size(y_node);
% size(z_node);
% figure
% plot(colvec(x_node),colvec(y_node),'k+');

% pressure
%pnode = mphinterp(model,{'p'},'coord',coords,'t',[PST.t1 PST.t2]); % pressure in Pa

% plot values at observation points
pts_coords=zeros(3,numel(DST.x));
%needed if not UTM already in Comsol so all models but the one with topo!!
%pts_coords(1,:)=DST.x-PST.xcen; % x coordinates in Comsol frame in meters
%pts_coords(2,:)=DST.y-PST.ycen; % y coordinates in Comsol frame in meters
%pts_coords(3,:)=zeros(size(DST.z)); % top of model is at z=0m

% when using the block with topography the UTM coordinates are already in
% Comsol

% pts_coords(1,:)=DST.x-PST.xcen; % x coordinates in Comsol frame in meters
% pts_coords(2,:)=DST.y-PST.ycen; % y coordinates in Comsol frame in meters
% xcen = PST.p0(get_parameter_index('CS_Easting_in_m_________________',PST.names))
% ycen = PST.p0(get_parameter_index('CS_Northing_in_m________________',PST.names))
% F#  95 CS_Origin_Easting_in_m__________        0.0        NaN        0.0        NaN     NaN        0.0
% F#  96 CS_Origin_Northing_in_m_________        0.0        NaN        0.0        NaN     NaN        0.0

ix = get_parameter_index('CS_Origin_Easting_in_m__________',PST.names);
if ix > 0
    xcen = PST.p0(ix);
else
    error('Cannot find CS_Origin_Easting_in_m__________ in parameter list\n');
end
iy = get_parameter_index('CS_Origin_Northing_in_m_________',PST.names);
if iy > 0
    ycen = PST.p0(iy);
else
    error('Cannot find CS_Origin_Northing_in_m_________ in parameter list\n');
end

if abs(xcen) < 1.0
    error(sprintf('xcen is small.'));
end
if abs(ycen) < 1.0
    error(sprintf('ycen is small.'));
end

pts_coords(1,:)=DST.x-xcen;
pts_coords(2,:)=DST.y-ycen;
%pts_coords(3,:)=DST.z; % top of model is the topographic surface (DEM values)
pts_coords(3,:)=zeros(size(DST.z)); % top of model is at z=0m
Emi=max(max(pts_coords(1,:)));
Ema=min(min(pts_coords(1,:)));
% figure;
% hist(pts_coords(1,:))

if numel(solution_epochs) > 1
    % for time-dependent viscoelastic solutions, we need interpolation in time
    % get absolute epochs in years
    t0 = PST.p0(get_parameter_index('Reference_Epoch_in_years________',PST.names));
    t1 = PST.p0(get_parameter_index('time_fn_@_epoch_001_in_years____',PST.names));
    t2 = PST.p0(get_parameter_index('time_fn_@_epoch_002_in_years____',PST.names));
    
    % convert to relative time in seconds 
    % with respect to initial conditions
    secperyr = 365.25 * 3600 * 24;
    dt1 = (t1 - t0) * secperyr;
    dt2 = (t2 - t0) * secperyr;
    % displacement at master epoch at observation points
    upts1 = colvec(mphinterp(model,{'u'},'coord',pts_coords,'t',dt1)); % easting component
    vpts1 = colvec(mphinterp(model,{'v'},'coord',pts_coords,'t',dt1)); % northing component
    wpts1 = colvec(mphinterp(model,{'w'},'coord',pts_coords,'t',dt1)); % vertical component
    % displacement at slave epoch at observation points
    upts2 = colvec(mphinterp(model,{'u'},'coord',pts_coords,'t',dt2)); % easting component
    vpts2 = colvec(mphinterp(model,{'v'},'coord',pts_coords,'t',dt2)); % northing component
    wpts2 = colvec(mphinterp(model,{'w'},'coord',pts_coords,'t',dt2)); % vertical component
    
    %differential displacements
    dupts=upts2-upts1;
    dvpts=vpts2-vpts1;
    dwpts=wpts2-wpts1;    
else
    % differential displacement at observation points - for elastic case
%     dupts = colvec(mphinterp(model,{'u'},'coord',pts_coords,'ext',1.0)); % easting component
%     dvpts = colvec(mphinterp(model,{'v'},'coord',pts_coords,'ext',1.0)); % northing component
%     dwpts = colvec(mphinterp(model,{'w'},'coord',pts_coords,'ext',1.0)); % vertical component
    dupts = colvec(mphinterp(model,{'u'},'coord',pts_coords)); % easting component
    dvpts = colvec(mphinterp(model,{'v'},'coord',pts_coords)); % northing component
    dwpts = colvec(mphinterp(model,{'w'},'coord',pts_coords)); % vertical component    
end

if verbose == 1
    fprintf(1,'Extrema in Differential Eastward   %12.4e %12.4e\n',nanmin(nanmin(dupts)),nanmax(nanmax(dupts)));
    fprintf(1,'Extrema in Differential Northward  %12.4e %12.4e\n',nanmin(nanmin(dvpts)),nanmax(nanmax(dvpts)));
    fprintf(1,'Extrema in Differential Upward     %12.4e %12.4e\n',nanmin(nanmin(dwpts)),nanmax(nanmax(dwpts)));
    
    hist(colvec(hypot(dupts,dvpts)));
    xlabel('Displacement[m]');
    ylabel('Number of points');
    
    size(DST.x)
    size(dupts)
    size(DST.y)
    size(dvpts)
    
    figure;
    quiver(colvec(DST.x/1e3),colvec(DST.y/1e3),colvec(dupts),colvec(dvpts));
end



return

end
