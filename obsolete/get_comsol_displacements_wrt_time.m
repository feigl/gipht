function [solution_epochs,xpts,ypts,zpts,dupts,dvpts,dwpts] = get_comsol_displacements_wrt_time(mphfile,pts_coords,tepochs)
%function [solution_epochs,xpts,ypts,zpts,dupts,dvpts,dwpts] = get_comsol_displacements(PST,DST,TST)
% Kurt 2014-08-08

verbose = 0;
%verbose = 1;

%[pnames32, pnames, values, dims, descrs] = get_comsol_parameters(PST.datafilename)

% Comsol evaluation solution_epochs in seconds
%[model, solution_epochs] = execute_comsol(PST.datafilename)

model = mphload(mphfile);

info1 = mphsolinfo(model);
info2 = mphsolutioninfo(model);

if verbose == 1
    info1
    info2
end

% Comsol evaluation solution_epochs in seconds
solution_epochs = info1.solvals;


%obtain nodal coordinates first:
nodestruct = mphxmeshinfo(model,'soltag','sol1');
coords = nodestruct.nodes.coords;


% pts_coords = zeros(3,1);
% pts_coords(1,1)=0.;
% pts_coords(2,1)=0.;
% pts_coords(3,1)=0.; 
[nr,nc] = size(pts_coords);
if nr ~= 3
    error
end
xpts = pts_coords(1,:);
ypts = pts_coords(2,:);
zpts = pts_coords(3,:);

if numel(solution_epochs) > 1
    % for time-dependent viscoelastic solutions, we need interpolation in time
    % get absolute epochs in years
     
    % convert to relative time in seconds 
    % with respect to initial conditions
    secperyr = 365.25 * 3600 * 24;

    
    dt1 = tepochs * secperyr
    % displacement at master epoch at observation points
    upts1 = mphinterp(model,'u','coord',pts_coords,'t',dt1); % easting component
    vpts1 = mphinterp(model,'v','coord',pts_coords,'t',dt1); % northing component
    wpts1 = mphinterp(model,'w','coord',pts_coords,'t',dt1); % vertical component
    % displacement at slave epoch at observation points
    upts2 = zeros(size(upts1)); % easting component
    vpts2 = zeros(size(vpts1)); % northing component
    wpts2 = zeros(size(wpts1)); % vertical component
    
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



return

end
