function DST = build_dst(fitfun,xyzm,tepochs,bpest,dops,DD,unitv...
    ,xd,yd,ippix1,mpercy,idatatype...
    ,dx,dy,dz,orbvm,orbvs,alon,alat...
    ,qii1,qii2,qjj1,qjj2...
    ,phasig)
% write all the data needed for a model to a data structure called DST

if nargin < 19
    error('Not enough arguments');
end

% get dimensions
[np,me] = size(DD);
[ndum,ndata] = size(xyzm);


% reserve space for the model
phamod = 0.;

if exist('orbvm','var') ~= 1
    warning('Missing orbvm');
    orbvm = zeros(6,ndata);
end
if exist('orbvs','var') ~= 1
    warning('Missing orbvs');
    orbvs = zeros(6,ndata);
end
if exist('alon','var') ~= 1
    warning('Missing alon');
    alon = nan(ndata,1);
end
if exist('alat','var') ~= 1
    warning('Missing alat');
    alat = nan(ndata,1);
end

% remove means
for i=1:3
    xyzm0(i,:) = xyzm(i,:)-mean(xyzm(i,:));
end

kmasts = zeros(ndata,1);
kslavs = zeros(ndata,1);
for k=1:np
    i1 = ippix1(k);
    if k < np
        %i2 = ippix1(k+1);
        i2 = ippix1(k+1)-1;
    else
        i2 = ndata;
    end
    npixinpair=i2-i1+1;
    
    kmast = find(DD(k,:) == -1);
    kslav = find(DD(k,:) == +1);
    
    kmasts(i1:i2) = kmast;
    kslavs(i1:i2) = kslav;
    kindex(i1:i2) = k;    
end

DST.i            = [1:ndata]';
DST.idatatype    = idatatype*ones(ndata,1);
DST.k            = colvec(kindex);
DST.kmast        = colvec(kmasts);
DST.kslav        = colvec(kslavs);
DST.phaobs       = colvec(xd);
DST.phamod       = zeros(ndata,1);
DST.tmast        = colvec(tepochs(kmasts));
DST.tslav        = colvec(tepochs(kslavs));
DST.x            = colvec(xyzm(1,:));
DST.y            = colvec(xyzm(2,:));
DST.z            = colvec(xyzm(3,:));
DST.uvx          = colvec(unitv(1,:));
DST.uvy          = colvec(unitv(2,:));
DST.uvz          = colvec(unitv(3,:));
% 2012-JUN-25 vector values for fringe interval (mpercy = meters per cycle)
if numel(mpercy) == ndata
    DST.mpercy       = colvec(mpercy);
elseif numel(mpercy) == 1
    DST.mpercy       = mpercy*ones(ndata,1);
else
    warning('Mismatch on mpercy\n');
    DST.mpercy       = nanmean(mpercy)*ones(ndata,1);
end
DST.bmast        = colvec(bpest(kmasts));
DST.bslav        = colvec(bpest(kslavs));
DST.dmast        = colvec(dops(kmasts));
DST.dslav        = colvec(dops(kslavs));
%New items 2010-OCT-05
DST.x0           = colvec(xyzm0(1,:));         % de-meaned coordinates
DST.y0           = colvec(xyzm0(2,:));         % de-meaned coordinates
DST.z0           = colvec(xyzm0(3,:));         % de-meaned coordinates
%These MUST BE CONSTANT. Their values come from the DEM descriptor
DST.dx           = colvec(dx*ones(ndata,1));   % easting  pixel dimension in meters
DST.dy           = colvec(dy*ones(ndata,1));   % northing pixel dimesnion in meters
DST.dz           = colvec(dz);                 % change of topography in eastward direction in meters
% New items 2011-JUN-26
% DST.bradi        = colvec(orbvm(1,:)); % component of Baseline vector parallel to radius through satellite
% DST.bhori        = colvec(orbvm(2,:)); % component of Baseline vector parallel to satellite velocity vector
% DST.bperp        = colvec(orbvm(3,:)); % component of Baseline vector perpendicular to line of sight
% DST.bpara        = colvec(orbvm(4,:)); % component of Baseline vector parallel to line of sight from sat to target
% DST.theta        = colvec(orbvm(5,:)); % look angle in radians (NOT the same as incidence)
% DST.ndist        = colvec(vnear(1,:)); % range from satellite at closest approach to target in meters
% DST.tnear        = colvec(vnear(2,:)); % time of closest approach to target in seconds of MJD
% DST.incid        = colvec(vnear(3,:)); % incidence angle at target in radians from vertical         
% DST.alonr          = colvec(alon); % geodetic longitude on WGS84 ellipsoid in radians    
% DST.alatr          = colvec(alat); % geodetic latitude  on WGS84 ellipsoid in radians  
DST.alond          = colvec(alon); % geodetic longitude on WGS84 ellipsoid in degrees 
DST.alatd          = colvec(alat); % geodetic latitude  on WGS84 ellipsoid in degrees  
DST.orbm1          = colvec(orbvm(1,:)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
DST.orbm2          = colvec(orbvm(2,:)); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
DST.orbm3          = colvec(orbvm(3,:)); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
DST.orbm4          = colvec(orbvm(4,:)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
DST.orbm5          = colvec(orbvm(5,:)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
DST.orbm6          = colvec(orbvm(6,:)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
DST.orbs1          = colvec(orbvs(1,:)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
DST.orbs2          = colvec(orbvs(2,:)); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
DST.orbs3          = colvec(orbvs(3,:)); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
DST.orbs4          = colvec(orbvs(4,:)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
DST.orbs5          = colvec(orbvs(5,:)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
DST.orbs6          = colvec(orbvs(6,:)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
% New items 2012-JUN-25
DST.mx           = nan(ndata,1);   % modeled displacement vector in meters: East   component
DST.my           = nan(ndata,1);   % modeled displacement vector in meters: North  component
DST.mz           = nan(ndata,1);   % modeled displacement vector in meters: Upward component
% New items 2012-OCT-04
DST.qii1         = colvec(double(qii1));   % lower limit of row index for quadtree patch
DST.qii2         = colvec(double(qii2));   % upper limit of row index for quadtree patch
DST.qjj1         = colvec(double(qjj1));   % lower limit of col index for quadtree patch
DST.qjj2         = colvec(double(qjj2));   % upper limit of col index for quadtree patch
% New items 2014-JAN-08
DST.phasig        = colvec(double(phasig)); % 1-sigma uncertainty of a single phase measurement (same units as phaobs)

 
return

