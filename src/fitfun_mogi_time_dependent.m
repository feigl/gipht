function Lmod = fitfun_mogi_time_dependent(PST,DST,TST)
% PST parameter structure
% DST data structure
% fitting function
%persistent kount

narginchk(3,3);

% %fprintf(1,'Entering %s\n',mfilename);
% mparams = numel(PST.p1);
% for i=1:mparams
% fprintf(1,'%10.4g ',PST.p1(i));
% end
% fprintf(1,'\n');

% times when data occur in seconds elapsed from start of model
nEpochs = numel(DST.t);

% number of data
ndata = numel(DST.obs);

% find unique locations by distance from orgin
rdist=sqrt(DST.x.^2 + DST.y.^2 + DST.z.^2);
[runique,iunique] = unique(rdist);
nObsPoints = numel(runique);

nObsPoints = numel(runique);
XYobspts=zeros(2,nObsPoints);

XYobspts(1,:) = reshape(DST.x(iunique),1,nObsPoints);
XYobspts(2,:) = reshape(DST.y(iunique),1,nObsPoints);


%tEpochs = reshape(DST.t - PST.t0,1,nEpochs);
%tEpochs = reshape(DST.t,1,nEpochs);
Lmod = nan(ndata,1);

%% calculate modeled values at locations of data points
volgeom(1) = PST.p1(1);          % E coordinate of source in m
volgeom(2) = PST.p1(2);          % N coordinate of source in m
volgeom(3) = PST.p1(3);          % Depth of source in m
nu         = PST.p1(5);          % Poisson ratio, dimensionless = 0.25;

Umod = nan(nEpochs,nObsPoints);
Vmod = nan(nEpochs,nObsPoints);
Wmod = nan(nEpochs,nObsPoints);

for i=1:nEpochs
    volgeom(4) = DST.t(i)*PST.p1(4);          % Volume change of source in m^3  = 1.382e+08;  
    umogi = mogi(volgeom, XYobspts, nu);
    Umod(i,:) = umogi(1,:);
    Vmod(i,:) = umogi(2,:);
    Wmod(i,:) = umogi(3,:);
end


% %% extract modeled displacements by interpolating coordinates and time
% Umod = mphinterp(model,{'u2'},'coord',XYZobpts,'t',tEpochs,'dataset',PST.dset); % easting component
% Vmod = mphinterp(model,{'v2'},'coord',XYZobpts,'t',tEpochs,'dataset',PST.dset); % northing component
% Wmod = mphinterp(model,{'w2'},'coord',XYZobpts,'t',tEpochs,'dataset',PST.dset); % vertical component


[nur,nuc] = size(Umod);
if nur ~= nEpochs || nuc ~= nObsPoints 
    nur
    nuc
    nEpochs
    nObsPoints
    warning('arrays are the wrong size');
end

% make displacement relative to first observed epoch 
for i=1:nObsPoints
    Umod(:,i) = Umod(:,i) - Umod(1,i);
    Vmod(:,i) = Vmod(:,i) - Vmod(1,i);
    Wmod(:,i) = Wmod(:,i) - Wmod(1,i);
end

% compute displacement along line of sight
Lmod = Umod.*DST.pointingX + Vmod.*DST.pointingY + Wmod.*DST.pointingZ;

return
end

