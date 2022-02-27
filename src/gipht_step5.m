% GIPHT_step5: make plots and images for entire sub-region
% 2009-DEC-03
% 2010-MAR-22 change totcost[00,0,1] to cost[00,0,1] to save calls to
% fitfun
% 2010-NOV-14 Use DST,PST,TST structures
% 2012-MAY-23
% 2012-JUN-27
% 2012-JUL-10 helene
% 2012-OCT-04 Add min and max values to PST file

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

clear vars;
if fexist('gipht.mat') == 1
    load('gipht.mat');
end

% initialize
kq = 0; % pointer to derived parameters

%load;
%load('qsave.mat');
fidtxtout = fopen(fnparout,'a');
% iq1 = iq;
% iq2 = iq;

% How handle values for multi-panel plots
%figopt % xx1 propagate nulls from quadtree, paint missing data black
%figopt % x1x re-calculate modeled values at all pixel locations
%figopt % x0x interpolate sample to find modeled values at all pixel locations
%figopt % 1xx request grids and profiles of vector components of displacement


%% Restore exact model instead of Taylor
% if ianneal == 4
%     fitfun = PST.fitfun;
% end
%
% %if  numel(strfind(fitfun,'funpressure4')) > 0
% pindex=get_parameter_index('annual',pnames);
% if numel(pindex) > 0
%     [tyr,ft,fts] = get_annual_rates(p1,psig,pnames,me,fitfun);
%     % remove points with zero uncertainties
%     iok = find(abs(fts)>0);
%     tyr=tyr(iok);
%     ft =ft(iok);
%     fts=fts(iok);
%
%     nf=nf+1;h(nf)=figure;hold on;
%     errorbar(tyr,ft,fts,'ko-');
%     xlabel('Year beginning');
%     ylabel('Modeled Annual Rate');
%     axis xy;
%     feval(printfun,sprintf('%s_Annual',runname));
% end

% put all models on one plot
%panelrows = ceil(sqrt(np));
%panelcols = floor(sqrt(np))
%panelcols = ceil(np/panelrows);
panelcols = 4;
panelrows = ceil(np/panelcols);
nfull = panelrows * panelcols;
climits = [-0.5 0.5];

% make a mesh of for entire sub-region
% ymax = max(yax(isub));
% ymin = min(yax(isub));
% xmax = max(xax(jsub));
% xmin = min(xax(jsub));
xmin = xsubmin;
xmax = xsubmax;
ymin = ysubmin;
ymax = ysubmax;
wesn = [xmin xmax ymin ymax];

disp 'Making meshes'
% Easting and Northing in meters
%[xmesh, ymesh]  = meshgrid(1e3*xax(jsub),1e3*yax(isub));
% mcells = round((xmax-xmin)/interpcell)
% ncells = round((ymax-ymin)/interpcell)
% mcells = numel(xax(jsub));
% ncells = numel(yax(isub));
mcells = numel(grdx);
ncells = numel(grdy);
if dl > 0
    [xmesh, ymesh]  = meshgrid(linspace(xmin,xmax,mcells),linspace(ymin,ymax,ncells)); % 1-km mesh
else
    [xmesh, ymesh]  = meshgrid(linspace(xmin,xmax,mcells),linspace(ymax,ymin,ncells)); % 1-km mesh
end

[nlmesh,ncmesh] = size(xmesh);
% bx = colvec(xmesh);
% by = colvec(ymesh);
bx = rowvec(xmesh);
by = rowvec(ymesh);
%bz = zeros(size(bx));
%bz = reshape(zframe(isub,jsub),1,nlmesh*ncmesh);
%bz = reshape(zframe,1,nlmesh*ncmesh);
bz = demz;
%
% % fine grid for all pixels in subregion
%[xmesh3,ymesh3] =  meshgrid(xax(jsub),yax(isub)); % all pixels in subregion
[xmesh3,ymesh3] =  meshgrid(grdx,grdy); % all pixels in subregion
[nlmesh3,ncmesh3] = size(xmesh3);
cy = colvec(ymesh3);
cx = colvec(xmesh3);

% Xmesh=reshape(bx,nlmesh,ncmesh);
% Ymesh=reshape(by,nlmesh,ncmesh);
% Zmesh=reshape(bz,nlmesh,ncmesh);
disp 'Done making meshes'

%% initialize unit vector
unitv1 = nan(3,ndata);

%% coordinates for profile
xprof = (xmin+xmax)/2.;
yprof = (ymin+ymax)/2.;
iprof = find(abs(demy-yprof) < dy);iprof=iprof(1);
jprof = find(abs(demx-xprof) < dx);jprof=jprof(1);


%% dimension arrays to store images as panels in a square mosaic
if np > 1 && np < 36
    oims = zeros(np,nlmesh3,ncmesh3);
    mims = zeros(np,nlmesh3,ncmesh3);
    rims = zeros(np,nlmesh3,ncmesh3);
    cims = zeros(np,nlmesh3,ncmesh3);
    tls  = cell(np,1);
end

% rename arrays 20150520
mCdl0 = mdl0; % modeled values at sample locations calculated from initial estimate
mCdl1 = mdl1; % modeled values at sample locations calculated from final estimate

disp('Reading phase files....')
phao_samples = DST.phaobs;
clear phao
% 20120926 cmapblackzero;
ctab = cmapblackzero(1); % black at zero at bottom of color bar

i2=1;j2=1;
for i = 1:np
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(1, 'making figures for pair %03d\n',i);
    
    clear imA imB imC imD imE imF imG imH;
    clear uns ucs urs mds;
    
    % select row corresponding to this pair
    %DD1 = DD(i,:);
    
    jmast =  abs(imast(i));
    jslav =  abs(islav(i));
    
    kmast = find(DD(i,:) == -1);
    kslav = find(DD(i,:) == +1);
    
    % select row corresponding to this pair
    DD1 = DD(i,:);
    
    %npixinpair = nlmesh3*ncmesh3;
    
    % New scheme using structures 2010-11-08
    kpairs = find(DST.k == i);
    j1 = DST.i(kpairs(1));
    j2 = DST.i(kpairs(end));
    % 20160524
    %time_span = DD1 * colvec(tepochs); % in years
    time_span = years(tepochs(kslav)-tepochs(kmast));
    
    fprintf(1,'Starting, stopping indices j1,j2 %d %d\n',j1,j2);
    
    %% read observed phase values at all pixels for this pair
    fn0 = pfnames{i};
    % Kurt 20170724 dimensions of subregion must be dimensions of input
    % grid file. In other words, cut the grid files down to size before
    % running GIPhT.
    INFO = grdinfo3(fn0);
    nrsub = INFO.ny
    ncsub = INFO.nx
    %[tmpx,tmpy,tmpz] = grdread3(fn0); % GMT grid file 
    [grdx,grdy,tmpz] = grdread3(fn0); % GMT grid file 
    phao = double(tmpz)*FACTIN;
    ndata = numel(phao);
    
    kmasts = kmast*ones(ndata,1);
    kslavs = kslav*ones(ndata,1);

    
    if bitget(figopt,2) == 1
        % recalculate model at each pixel individually
        i1=1;
        i2=nlmesh*ncmesh;
        ndata = nlmesh*ncmesh;
        
        ippix1(1)=i1;
        ippix1(2)=i2;
        
        clear xyzm1; clear unitv1;
        xyzm1(1,i1:i2) = rowvec(bx);
        xyzm1(2,i1:i2) = rowvec(by);
        xyzm1(3,i1:i2) = rowvec(bz);
        
        % zero some values
        %         dx1(i1:i2) = zeros(1,ndata); %
        %         dy1(i1:i2) = zeros(1,ndata); %
        %         dz1(i1:i2) = zeros(1,ndata); % topographic relief set to zero CAREFUL
        % zero some values
        
        %% Step sizes equal differences in coordinates, assumed constant
        dx1 = grdx(2)-grdx(1);
        dy1 = grdy(2)-grdy(1);
        dz1 = 0.; % topographic relief set to zero CAREFUL
        
        fprintf(1,'In routine %s matrix of coordinates bx    has dimensions %d\n',mfilename,numel(bx));
        fprintf(1,'In routine %s matrix of coordinates xyzm1 has dimensions %d %d\n',mfilename,size(xyzm1));
    else
        % interpolate values from sample - works fine for most
        % applications with sufficient number of sample
        DST2 = DST;
        % indices into images of subregion
        i1 = ippix1(i);
        if i < np
            i2 = ippix1(i+1)-1;
        else
            i2 = ndata;
        end
        
        ippix1(1)=i1;
        ippix1(2)=i2;
    end
    
    if numel(orbfile) > 0
        %         % individual orbit file
        %         fprintf(1,'Calculating look vectors for pair %d from orbit file %s\n',i,ofnames1{i});
        %         clear unitv orbvm orbvs;
        %
        %         % % coordinates of corners in meters for corners
        %         ux(1) = xmin;    uy(1) = ymin; % SW
        %         ux(2) = xmax;    uy(2) = ymin; % SE
        %         ux(3) = xmax;    uy(3) = ymax; % NE
        %         ux(4) = xmin;    uy(4) = ymax; % NW
        %         ux(5) = xcenter; uy(5) = ycenter;     % center
        %         for ii=1:5
        %             xyzm1(1,ii) = ux(ii);
        %             xyzm1(2,ii) = uy(ii);
        %             xyzm1(3,i)  = 0;
        %             phao1(ii)   = 0.0;
        %             dx1(ii)     = 0;
        %             dy1(ii)     = 0;
        %             dz1(ii)     = 0;
        %         end
        %
        %         % individual orbit file
        %         for ii=1:5
        %             % approximate latitude and longitude from inverse UTM
        %             if isgeo == 2
        %                 [blat,blon]=utm2deg(ux(ii)-x0+xutm0,uy(ii)-y0+yutm0,utmzone0);
        %             else
        %                 if ii==1 || ii==4
        %                     blon = min(alon);
        %                 elseif ii==2 || ii==3
        %                     blon = max(alon);
        %                 else
        %                     blon = mean([min(alon), max(alon)]);
        %                 end
        %                 if ii==1 || ii==4
        %                     blat = min(alat);
        %                 elseif ii==2 || ii==3
        %                     blat = max(alat);
        %                 else
        %                     blat = mean([min(alat), max(alat)]);
        %                 end
        %             end
        %             %             % calculate unit vector that varies across scene
        %             %             [uv1(ii),uv2(ii),uv3(ii)] = lookvector(blon,blat,0.,ofnames1{i});
        %             fprintf(1,'Calculating orbital information at lon = %#12.4f lat %#12.4f for pair %d from orbit files %s and %s\n'...
        %                 ,blon,blat,i,ofnames1{i},ofnames2{i});
        %
        % %             % calculate unit vector that varies across scene
        % %             [uv1(ii),uv2(ii),uv3(ii),near_mjd,near_sec(ii),near_dist(ii),near_vel, incid(ii)] ...
        % %                 = lookvector(blon,blat,0.,ofnames1{i});
        % %
        % %             % calculate components of baseline compnents in meters
        % %             [bradi(ii), bhori(ii), bperp(ii), bpara(ii), theta(ii)] ...
        % %                 = orbdiff(blon,blat,0.,ofnames1{i},ofnames2{i});
        %
        %              [UT, DT, NT, RT, HT, AT, VT, near_mjdT, near_secT] = orbvectors2(loncenter,latcenter,0.0 ...
        %                 ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
        %             uv1(ii) = UT(1);
        %             uv2(ii) = UT(2);
        %             uv3(ii) = UT(3);
        %         end
        
        %         unitv1(1,i1:i2) = griddata(ux,uy,uv1,bx,by,'v4'); % east  components
        %         unitv1(2,i1:i2) = griddata(ux,uy,uv2,bx,by,'v4'); % north components
        %         unitv1(3,i1:i2) = griddata(ux,uy,uv3,bx,by,'v4'); % up components
        %
        % %         vnear(1,i1:i2) = griddata(ux,uy,near_dist,bx,by,'v4'); % range from satellite at closest approach to target in meters
        % %         vnear(2,i1:i2) = griddata(ux,uy,near_sec ,bx,by,'v4');  % time of closest approach to target in seconds of MJD
        % %         vnear(3,i1:i2) = griddata(ux,uy,incid    ,bx,by,'v4');     % incidence angle at target in radians from vertical
        % %         orbvs(1,i1:i2) = griddata(ux,uy,bradi    ,bx,by,'v4'); % component of Baseline vector parallel to radius through satellite
        % %         orbvs(2,i1:i2) = griddata(ux,uy,bhori    ,bx,by,'v4'); % component of Baseline vector parallel to satellite velocity vector
        % %         orbvs(3,i1:i2) = griddata(ux,uy,bperp    ,bx,by,'v4'); % component of Baseline vector perpendicular to line of sight
        % %         orbvs(4,i1:i2) = griddata(ux,uy,bpara    ,bx,by,'v4'); % component of Baseline vector parallel to line of sight from sat to target
        % %         orbvs(5,i1:i2) = griddata(ux,uy,theta    ,bx,by,'v4'); % look angle in radians
        %         orbvm1(1:6,i1:i2) = 0;
        %         orbvs1(1:6,i1:i2) = 0;
        %         te0=tic; unitv1(1,i1:i2) = griddata(ux,uy,unitv(1,:),bx,by,'v4'); fprintf(1,'%.1f seconds\n',toc(te0));
        %         te0=tic; unitv1(2,i1:i2) = griddata(ux,uy,unitv(2,:),bx,by,'v4'); fprintf(1,'%.1f seconds\n',toc(te0));
        %         te0=tic; unitv1(3,i1:i2) = griddata(ux,uy,unitv(3,:),bx,by,'v4'); fprintf(1,'%.1f seconds\n',toc(te0));
        
        
        
        %         for jj = j1:j2
        %             fprintf(1,'%#16.1f %#16.1f\n',DST.x(jj),DST.y(jj));
        %         end
        
        if bitget(figopt,2) == 1
            % 20130629 use new griddata2 instead of old griddata
            fprintf(1,'Interpolating 3 components of unit vectors....\n');
            unitv1 = nan(3,numel(bx));
            unitv1(1,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.uvx(j1:j2),bx,by,'linear');
            unitv1(2,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.uvy(j1:j2),bx,by,'linear');
            unitv1(3,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.uvz(j1:j2),bx,by,'linear');
            fprintf(1,'Interpolating 6 orbit vectors for master....\n');
            orbvm1 = zeros(6,numel(bx));
            orbvm1(1,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm1(j1:j2),bx,by,'linear');
            orbvm1(2,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm2(j1:j2),bx,by,'linear');
            orbvm1(3,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm3(j1:j2),bx,by,'linear');
            orbvm1(4,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm4(j1:j2),bx,by,'linear');
            orbvm1(5,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm5(j1:j2),bx,by,'linear');
            orbvm1(6,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbm6(j1:j2),bx,by,'linear');
            fprintf(1,'Interpolating 6 orbit vectors for slave...\n');
            orbvs1 = nan(6,numel(bx));
            orbvs1(1,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs1(j1:j2),bx,by,'linear');
            orbvs1(2,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs2(j1:j2),bx,by,'linear');
            orbvs1(3,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs3(j1:j2),bx,by,'linear');
            orbvs1(4,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs4(j1:j2),bx,by,'linear');
            orbvs1(5,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs5(j1:j2),bx,by,'linear');
            orbvs1(6,:) = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.orbs6(j1:j2),bx,by,'linear');
            fprintf(1,'Interpolating latitudes and longitudes of pixel centers.\n');
            % geodetic longitude in degrees
            clond1 = nan(1,numel(bx));
            clond1(1,:)  = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.alond(j1:j2),bx,by,'linear');
            % geodetic latitude in degrees
            clatd1 = nan(1,numel(bx));
            clatd1(1,:)  = griddata2(DST.x(j1:j2),DST.y(j1:j2),DST.alatd(j1:j2),bx,by,'linear');
        end
    else
        % constant across scene
        unitv1(1,i1:i2)=unitv0(1);
        unitv1(2,i1:i2)=unitv0(2);
        unitv1(3,i1:i2)=unitv0(3);
        orbvm1=zeros(6,ndata);
        orbvs1=zeros(6,ndata);
        fprintf(1,'Zeroing latitudes and longitudes of pixel centers.\n');
        clond1 = zeros(1,ndata);
        clatd1 = zeros(1,ndata);
    end
    
    % set up DST2 structure for interpolation
    if bitget(figopt,2) == 1
        % calculate modeled values for phase, NOT gradient
        %idatatype = 0;
        % calculate modeled values for whatever we are using
        idatatype = idatatype1;
        
        storage=[];
        clear DST2;
%         function DST = build_dst(fitfun,xyzm,tepochs,bpest,dops,DD,unitv...
%     ,xd,yd,ippix1,mpercys,idatatype...
%     ,dx,dy,dz,orbvm,orbvs,alon,alat...
%     ,qii1,qii2,qjj1,qjj2...
%     ,phasig,kindex,kmasts,kslavs)

        DST2 = build_dst(fitfun,xyzm1,tepochs,bpest,dops,DD1,unitv1...
            ,zeros(ndata,1),zeros(ndata,1),ippix1,mpercys(i),idatatype...
            ,dx1,dy1,dz1,orbvm1,orbvs1,clond1,clatd1...
            ,qii1,qii2,qjj1,qjj2...
            ,phasig,kindex,kmasts,kslavs);
        
        fprintf(1,'Number of finite values in DST2.z0 = %d\n',numel(find(isfinite(DST2.z0)==1)));
    
        %             % call fitting function first time to initialize
        %             % calculate modeled values for phase, NOT gradient
        for ii=1:numel(DST2.idatatype)
            DST2.idatatype(ii) = idatatype1;
        end
        
        % use exact fitting function
        fitfun = fitfun_exact;
        % temporary storage structure TST
        %[rng,TST] = feval(fitfun,DST2,PST);
        [rng,TST0] = feval(fitfun,DST2,PST0); %20150727
        % initial model
        fprintf(1,'Evaluating fitting function %s for initial estimate at %d locations....\n',fitfun,numel(DST2.x));
        %mCdl0 = feval(fitfun,DST2,PST,TST); % initl
        mCdl0 = feval(fitfun,DST2,PST0,TST0); % 20150727
        if isreal(mCdl0) ~= 1
            warning(sprintf('Found complex values in mCdl0. Max abs(imaginary) is %g\n',max(imag(mCdl0))));
            mCdl0 = real(mCdl0);
        end
        %             fprintf(1,'Setting Offset parameters in final estimate to zero ....\n');
        %             PST1.p1(22) = 0;
                % temporary storage structure TST
        %[rng,TST] = feval(fitfun,DST2,PST);
        [rng,TST1] = feval(fitfun,DST2,PST1); %20150727

        fprintf(1,'Evaluating fitting function %s for final   estimate at %d locations....\n',fitfun,numel(DST2.x));
        mCdl1 = feval(fitfun,DST2,PST1,TST1); % final
        %profile off
        if isreal(mCdl1) ~= 1
            warning(sprintf('Found complex values in mCdl0. Max abs(imaginary) is %g\n',max(imag(mCdl0))));
            mCdl1 = real(mCdl1);
        end
        
        % phase residuals for all pixels in this pair
        nf=nf+1;h(nf)=figure;
        res1 = rwrapm(DST.phaobs-DST.phamod);
%         % Test for von Miseness
%         [Sm,VMnessString1] = vonmisesness (colvec(res1))    % Mardia and Jupp
%         [U2,VMnessString2] = vonmisesness2(colvec(res1))    % Fisher
%         % make QQ plot
%         qqplotvonmises(colvec(res1)/2.0/pi);
%         feval(printfun,sprintf('%s_QQPLOT_%03d_von_mises',runname,i));
        
        % get vector values
        if bitget(figopt,3) == 1
            fprintf(1,'Evaluating fitting function %s for final   estimate at %d locations....\n',fitfun,numel(DST2.x));
           %KF20170724[rng0,DST2] = feval(fitfun,DST2,PST1,TST); % final
            [rng0,DST2] = feval(fitfun,DST2,PST1,TST1); % final
            %profile off
            if isreal(rng0) ~= 1
                warning(sprintf('Found complex values in rng0. Max abs(imaginary) is %g\n',max(imag(rng0))));
                rng0 = real(rng0);
            end
            fprintf(1,'Number of defined model vectors %d\n',numel(find(isfinite(DST2.mx))));
        end
        % make array of modeled values into an image with same dimensions as observed image.
        mdl0 = reshape(mCdl0,size(phao));
        mdl1 = reshape(mCdl1,size(phao));
    else
        % interpolate model values at all pixel locations in subregion
        
        fprintf(1,'Evaluating interpolation functions in subregion with %d lines and %d columns \n',nlmesh3,ncmesh3);
        ti0 = tic;
        % evaluate interpolating functions at pixel locations
        % resample onto pixel locations
        % %2011-NOV-30 Out of memory. Type HELP MEMORY for your options.
        % %2011-NOV-30               mdl0 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'v4');
        % %2011-NOV-30               mdl1 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'v4');
        %     mdl0 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'cubic');
        %     mdl1 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'cubic');
        %     mdl0 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'nearnat');
        %     mdl1 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'nearnat');
        %     mdl0 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'linear');
        %     mdl1 = griddata(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'linear');
        %             F0 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),'natural');  % makes zits
        %             F1 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),'natural');
        %          % 20120926 - use nearest neighbor interpolation because it returns values at the edges
        %         F0 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),'nearnat');  % looks strange
        %         F1 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),'nearnat');
        %         mdl0=F0(xmesh3,ymesh3); %
        %         mdl1=F1(xmesh3,ymesh3); %
        %
        % %         F0 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),'linear');
        % %         F1 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),'linear');
        %         F0 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),'natural');
        %         F1 = TriScatteredInterp(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),'natural');
        %         mdl0=F0(xmesh3,ymesh3); %
        %         mdl1=F1(xmesh3,ymesh3); %
        %         mdl0 = griddata2(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'linear');
        %         mdl1 = griddata2(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'linear');
        mdl0 = griddata2(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl0(i1:i2)),xmesh3,ymesh3,'nearnat');
        mdl1 = griddata2(colvec(DST.x(i1:i2)),colvec(DST.y(i1:i2)),colvec(mCdl1(i1:i2)),xmesh3,ymesh3,'nearnat');
        
        %         mdl0=F0(cx,cy); % row vector
        %         mdl1=F1(cx,cy);  % row vector
        %             Fmx = TriScatteredInterp(colvec(DST2.x(i1:i2)),colvec(DST2.y(i1:i2)),colvec(DST2.mx(i1:i2)),'linear');
        %             Fmy = TriScatteredInterp(colvec(DST2.x(i1:i2)),colvec(DST2.y(i1:i2)),colvec(DST2.my(i1:i2)),'linear');
        %             Fmz = TriScatteredInterp(colvec(DST2.x(i1:i2)),colvec(DST2.y(i1:i2)),colvec(DST2.mz(i1:i2)),'linear');
        %             umx = Fmx(xmesh3,ymesh3); % modeled displacement vector in m [east    component]
        %             umy = Fmy(xmesh3,ymesh3); % modeled displacement vector in m [north   component]
        %             umz = Fmz(xmesh3,ymesh3); % modeled displacement vector in m [upward  component]
        
        %             disp 'mdl0'; size(mdl0)
        %             disp 'mdl1'; size(mdl1)
        fprintf(1,'Finished evaluating interpolation functions in %#10.4f seconds\n',toc(ti0));
        fprintf(1,'Interpolated grid of models have %d rows and %d columns\n',size(mdl0,1),size(mdl0,2));
        %     elseif all(size(mCdl0) == size(phao)) == 1 ...
        %             && all(size(mCdl1) == size(phao)) == 1 ...
        %             && bitget(figopt,2) == 1
    end
    
    % add mean of observed values
    if ismember(pselect,[7,9]) == 1
        %         mdl0 = mdl0 + rwrapm(mean_direction(phao) - mean_direction(mdl0));
        %         mdl1 = mdl1 + rwrapm(mean_direction(phao) - mean_direction(mdl1));
        %         mdl0 = mdl0 + nanmean(colvec(phao));
        %         mdl1 = mdl1 + nanmean(colvec(phao));
        mdl0 = mdl0 + mean_direction(phao);
        mdl1 = mdl1 + mean_direction(phao);
    end
    
    fprintf(1,'Extreme values of mdl0 %12.4e %12.4e \n',nanmin(nanmin(mdl0)),nanmax(nanmax(mdl0)));
    fprintf(1,'Extreme values of mdl1 %12.4e %12.4e \n',nanmin(nanmin(mdl1)),nanmax(nanmax(mdl1)));
    
    


    %% handle wrapping
    switch idatatype1
        case 0 %% observable is phase in radians
            %% residuals in radians
            res0 = rwrapm(phao-mdl0);
            res1 = rwrapm(phao-mdl1);
            
            % wrapped modeled values in radians
            wrm0 = rwrapm(mdl0);
            wrm1 = rwrapm(mdl1);
            %     add mean direction from obs
            %wrm0 = rwrapm(wrm0 + mean_direction(phao) - mean_direction(wrm0));
            %wrm0 = rwrapm(wrm0 + mean_direction(phao) - mean_direction(wrm0));
            %wrm0 = rwrapm(mdl0 + mean_direction(phao));
            % cannot estimate mean direction from gradients
            %if ismember(pselect,[5,7,9]) || pselect == 1
            %     if ismember(pselect,[1,5,7,9])
            %         wrm1 = rwrapm(wrm1 + mean_direction(phao) - mean_direction(wrm1));
            %     end
            
            fprintf(1,'Extreme values of wrm0 %12.4e %12.4e \n',nanmin(nanmin(wrm0)),nanmax(nanmax(wrm0)));
            fprintf(1,'Extreme values of wrm1 %12.4e %12.4e \n',nanmin(nanmin(wrm1)),nanmax(nanmax(wrm1)));
            
            % angular deviations for all pixels in subregion
            devs_all00 = rarcm(phao,zeros(size(phao)));
            devs_all0  = rarcm(phao,rwrapm(mdl0));
            devs_all1  = rarcm(phao,rwrapm(mdl1));
            
            % 2012-JUN-25 select OK points
            if pselect == 0
                iok = ones(size(devs_all00));
            else
                iok=find(devs_all00>0);
            end
            
             % 2010-MAR-22 OK to calculate scalar cost as average of vector costs
            totcost00 = nanmean(colvec(devs_all00(iok)));
            totcost0  = nanmean(colvec(devs_all0(iok)));
            totcost1  = nanmean(colvec(devs_all1(iok)));
            
            objfun1 = objfun;
        otherwise
            warning(sprintf('idatatype1 is %d. Skipping wrapping operation.\n',idatatype1));
            
            %% residuals in radians
            res0 = phao-mdl0;
            res1 = phao-mdl1;

            wrm0 = nan(size(mdl0));
            wrm1 = nan(size(mdl0));
            
            % for unwrapped data, angular deviations are equivalent to
            % residuals
            devs_all00 = phao;
            devs_all0  = res0;
            devs_all1  = res1;
                        
            % 2012-JUN-25 select OK points
            if pselect == 0
                iok = ones(size(devs_all00));
            else
                iok=find(isfinite(devs_all00) == 1);
            end

            % use RMS to evaluate cost
            totcost00 = nanrms(colvec(devs_all00(iok)));
            totcost0  = nanrms(colvec(devs_all0(iok)));
            totcost1  = nanrms(colvec(devs_all1(iok)));
            
            objfun1 = 'nanrms';
    end
    
 
    
    
    fprintf(1,'Total Cost (by %s) of null  model = %12.4f %s for %6d observations in sub-region\n',objfun1,totcost00,objlabel,numel(iok));
    fprintf(1,'Total Cost (by %s) of initl model = %12.4f %s for %6d observations in sub-region\n',objfun1,totcost0, objlabel,numel(iok));
    fprintf(1,'Total Cost (by %s) of final model = %12.4f %s for %6d observations in sub-region\n',objfun1,totcost1, objlabel,numel(iok));
    
    % Added 2008-JUL-17 Kurt
    xax2 = [xmin, xmax];
    yax2 = [ymax, ymin];
    
    
    % Co-seismic fault
    %%i=findparamindex(pnames,mparam,'Okada1_Length');
    %[Xcorners10,Ycorners10,Hcorners10,Ncorners] = disloc_to_seismo(p0(i:i+9));
    %[Xcorners11,Ycorners11,Hcorners11,Ncorners] = disloc_to_seismo(p1(i:i+9));
    % post-seismic fault
    %i=findparamindex(pnames,mparam,'Okada2_Length');
    %[Xcorners20,Ycorners20,Hcorners20,Ncorners] = disloc_to_seismo(p0(i:i+9));
    %[Xcorners21,Ycorners21,Hcorners21,Ncorners] = disloc_to_seismo(p1(i:i+9));
    
    %fprintf(1,'Coordinates of corners of Okada Fault patches\n')
    %for i=1:numel(Xcorners10)
    %   fprintf(1,'%20s %12.4f %12.4f\n',Ncorners{i},Xcorners11(i)/1000,Ycorners11(i)/1000);
    %end
    %for i=1:numel(Xcorners20)
    %   fprintf(1,'%20s %12.4f %12.4f\n',Ncorners{i},Xcorners21(i)/1000,Ycorners21(i)/1000);
    %end
    
    % draw corners of final model
    for ii=1:8
        mysyms{ii} = '';
        marksizes(ii) = 1;
    end
    for ii=[5 6 7]
        mysyms{ii} = 'ko-';
        marksizes(ii) = 2;
    end
    
    % observed phase for this pair
    phim =     reshape(double(phao)/DNPC,nlmesh3,ncmesh3); % convert from radians to cycles
    
    % unwrapped values in radians
    umd0 = mdl0 + res0;
    umd1 = mdl1 + res1;
    
    %% calculate unwrapped quantities in meters
    switch idatatype1
        case 0 %% observable is phase in radians
            if abs(mean(DST2.mpercy) - DST2.mpercy(1)) < 1.0e-6
                uns = DST2.mpercy(1)*reshape(     umd1/DNPC,nlmesh3,ncmesh3);  % unwrapped "observed" values
                mds = DST2.mpercy(1)*reshape(     mdl1/DNPC,nlmesh3,ncmesh3);  % unwrapped model
                urs = DST2.mpercy(1)*reshape(     res1/DNPC,nlmesh3,ncmesh3);  % unwrapped residuals
                ucs = DST2.mpercy(1)*reshape(devs_all1/DNPC,nlmesh3,ncmesh3);  % unwrapped costs (deviations)
            else
                error(sprintf('Inconsistent values of fringe spacing %20.8g %20.8g %20.8g\n'...
                    ,mean(DST2.mpercy),DST2.mpercy(1),mean(DST2.mpercy)-DST2.mpercy(1)));
                %         uns = 28.4*reshape(double(umd1)/DNPC,nlmesh3,ncmesh3);%uns=uns-mean(mean(uns));
                %         mds = 28.4*reshape(double(mdl1)/DNPC,nlmesh3,ncmesh3);%mds=mds-mean(mean(mds));
            end
        case -1 %% observable is dimensionless strain
            uns = reshape(     umd1,nlmesh3,ncmesh3);  % unwrapped "observed" values
            mds = reshape(     mdl1,nlmesh3,ncmesh3);  % unwrapped model
            urs = reshape(     res1,nlmesh3,ncmesh3);  % unwrapped residuals
            ucs = reshape(devs_all1,nlmesh3,ncmesh3);  % unwrapped costs (deviations)
        case 2 %% observable is range change in meters
            uns = reshape(     umd1,nlmesh3,ncmesh3);  % unwrapped "observed" values
            mds = reshape(     mdl1,nlmesh3,ncmesh3);  % unwrapped model
            urs = reshape(     res1,nlmesh3,ncmesh3);  % unwrapped residuals
            ucs = reshape(devs_all1,nlmesh3,ncmesh3);  % unwrapped costs (deviations)
        otherwise
            error(sprintf('Unknown datatype %d\n',idatatype1));
    end
    
    
    if bitget(figopt,3) == 1 && bitget(figopt,2) == 1
        %unwrapped components of displacement vectors in meters
        %function V  = range2vector(R,S,M)
        % given range, find vector
        %    inputs:
        %       r = scalar range change
        %       s = unit vector pointing from target to sensor
        
        % this one is smooth
        umr = mds; % 2012-JUN-26 WORKS!!
        figure;plot(colvec(DST2.x0),colvec(umr),'.k');title('umr');
        
        % get vector components, undefined below threshold
        %threshm = 0.5*cost1*mpercy; % depends on goodness of model
        threshm = 0.001; % constant, in meters
        %         umv = range2vector(rowvec(umr),[DST2.uvx';DST2.uvy';DST2.uvz'],[DST2.mx';DST2.my';DST2.mz'],threshm);
        % 2012-OCT-28 Much simpler
        umv = [DST2.mx';DST2.my';DST2.mz'];
        
        umr = reshape(umr     ,nrsub,ncsub); % unwrapped observed range  displacement in m
        umx = reshape(umv(1,:),nrsub,ncsub); % unwrapped observed vector displacement in m [east   component]
        umy = reshape(umv(2,:),nrsub,ncsub); % unwrapped observed vector displacement in m [north  component]
        umz = reshape(umv(3,:),nrsub,ncsub); % unwrapped observed vector displacement in m [upward component]
        fprintf(1,'Number of defined model vectors % d %d %d %d\n'...
            ,numel(find(isfinite(umv))) ...
            ,numel(find(isfinite(umx))) ...
            ,numel(find(isfinite(umy))) ...
            ,numel(find(isfinite(umz))));
        
        % unwrapped in meters
        omr = uns; % 2012-JUN-26 WORKS!!
        
        fprintf(1,'Before range2vector %d %d\n',size(omr));
        omv = range2vector(rowvec(omr),[DST2.uvx';DST2.uvy';DST2.uvz'],[DST2.mx';DST2.my';DST2.mz'],threshm);
        fprintf(1,'After  range2vector %d %d\n',size(omr));
        omr = reshape(omr,nrsub,ncsub);
        fprintf(1,'After  reshape      %d %d\n',size(omr));
%         figure;
%         plot(colvec(DST2.x0),colvec(omr),'.k');
%         title('omr');
        
        omx = reshape(omv(1,:),nrsub,ncsub); % unwrapped observed vector displacement in m [east   component]
        omy = reshape(omv(2,:),nrsub,ncsub); % unwrapped observed vector displacement in m [north  component]
        omz = reshape(omv(3,:),nrsub,ncsub); % unwrapped observed vector displacement in m [upward component]
        %         Fx = TriScatteredInterp(colvec(DST2.x),colvec(DST2.y),colvec(omv(1,:)),'linear');
        %         Fy = TriScatteredInterp(colvec(DST2.x),colvec(DST2.y),colvec(omv(2,:)),'linear');
        %         Fz = TriScatteredInterp(colvec(DST2.x),colvec(DST2.y),colvec(omv(3,:)),'linear');
        %         omx = Fx(xmesh3,ymesh3); % observed displacement vector in m [east    component]
        %         omy = Fy(xmesh3,ymesh3); % observed displacement vector in m [north   component]
        %         omz = Fz(xmesh3,ymesh3); % observed displacement vector in m [upward  component]
    else
        umr = NaN;
        umv = NaN;
        umx = NaN;
        umy = NaN;
        umz = NaN;
        omr = NaN;
        omv = NaN;
        omx = NaN;
        omy = NaN;
        omz = NaN;
    end
    
    %% title string for all plots
    titlestr = sprintf('Pair %3d epochs %3d %3d %s to %s Dt = %.1f yr '...
        ,i,iuniqorbs(kmast),iuniqorbs(kslav)...
        ,char(tepochs(kmast)),char(tepochs(kslav)),double(years(tepochs(kslav)-tepochs(kmast))));
      %  ,strrep(runname,'_','\_'));  
   %        ,nanmean(colvec(devs_all0))/DNPC,nanmean(colvec(devs_all1))/DNPC...
 
    
    %% build images
    imA = reshape(     double(phao)/DNPC ,nlmesh3,ncmesh3);tlA = 'Initial';
    if idatatype1 == 0
        imB = reshape(     double(wrm0)/DNPC ,nlmesh3,ncmesh3);tlB = 'Mod0';
    else
        imB = reshape(     double(mdl0)/DNPC ,nlmesh3,ncmesh3);tlB = 'Mod0';
    end
    imC = reshape(     double(res0)/DNPC ,nlmesh3,ncmesh3);tlC = 'Res0';
    imD = reshape(double(devs_all0)/DNPC ,nlmesh3,ncmesh3);tlD = 'Dev0';
    %if pselect == 3 || pselect == 5 || pselect == 7
    %if ismember(pselect,[3,5,7])
    if idatatype1 == -1
        %% TODO 20160814need to find quad-tree
        %imE = reshape(     double(qhao)/DNPC ,nlmesh3,ncmesh3);tlE = 'Final';
        imE = reshape(     double(phao)/DNPC ,nlmesh3,ncmesh3);tlE = 'Final';
    else
        imE = reshape(     double(phao)/DNPC ,nlmesh3,ncmesh3);tlE = 'Final';
    end
    if idatatype1 == 0       
        imF = reshape(     double(wrm1)/DNPC ,nlmesh3,ncmesh3);tlF = 'Mod1';
    else
        imF = reshape(     double(mdl1)/DNPC ,nlmesh3,ncmesh3);tlF = 'Mod1';        
    end
        
    imG = reshape(     double(res1)/DNPC ,nlmesh3,ncmesh3);tlG = 'Res1';
    imH = reshape(double(devs_all1)/DNPC ,nlmesh3,ncmesh3);tlH = 'Dev1';
    
    %   propagate zeros or missing data as black pixels
    % areas with elevation less than 1 m
    if nanstd(bz) > 1
        inull=find(abs(bz) < 1);
    else
        inull = [];
    end
    
%     
%     inull=union(inull,find(abs(imE) < 1/256.0));
%     inull=union(inull,find(isfinite(imE)==0)); %
%     inull=union(inull,find(abs(imA) < 1/256.0));
%     inull=union(inull,find(isfinite(imA)==0));
    
    % How handle values for multi-panel plots
    %figopt % xx1 propagate nulls from quadtree, paint missing data black
    %figopt % x1x calculate modeled values at all pixel locations
    %figopt % 1xx request grids and profiles of vector components of displacement
    
    if bitget(figopt,1) == 1
        %Pass nulls from observations to residuals and devs_all
        
        warning(sprintf('Replacing %d pixels with NaN\n',numel(inull)));
        % 20120926 leave models untouched
        imA(inull) = NaN;
        %imB(inull) = NaN;
        imC(inull) = NaN;
        imD(inull) = NaN;
        imE(inull) = NaN;
        %imF(inull) = NaN;
        imG(inull) = NaN;
        imH(inull) = NaN;
        % 2011-12-12 do unwrapped quantities too
        uns(inull) = NaN;
        urs(inull) = NaN;
        ucs(inull) = NaN;
        %mds(inull) = NaN;
    else
        % interpolate
        if idatatype1 == 0
            interpolation_method = 'nearest';extrapolation_method = 'none';
        else
            interpolation_method = 'natural';extrapolation_method = 'none';
        end

        warning(sprintf('Replacing %d pixels with Interpolated values using interpolation_method %s and extrapolation_method %s\n'...
            ,numel(inull)),interpolation_method,extrapolation_method);
        imA(inull) = NaN;
        %imA(inull) = fillingaps(imA,interpolation_method,extrapolation_method);
        %imB(inull) = NaN; 
        imC = fillingaps(imC,interpolation_method,extrapolation_method);
        imD = fillingaps(imD,interpolation_method,extrapolation_method);
        imE = fillingaps(imE,interpolation_method,extrapolation_method);
        imF = fillingaps(imF,interpolation_method,extrapolation_method);
        imG = fillingaps(imG,interpolation_method,extrapolation_method);
        imH = fillingaps(imH,interpolation_method,extrapolation_method);
        % do unwrapped quantities too
        uns = fillingaps(uns,interpolation_method,extrapolation_method);
        urs = fillingaps(urs,interpolation_method,extrapolation_method);
        ucs = fillingaps(ucs,interpolation_method,extrapolation_method);
        %mds(inull) = NaN;
    end
    
    % 20131106
    if dl > 0
        imA = flipud(imA);
        imB = flipud(imB);
        imC = flipud(imC);
        imD = flipud(imD);
        imE = flipud(imE);
        imF = flipud(imF);
        imG = flipud(imG);
        imH = flipud(imH);
    end
    
    % only if figopt == 111
    if bitget(figopt,1) == 1 && bitget(figopt,2) == 1 && bitget(figopt,3) == 1
        omx(inull) = NaN;
        omy(inull) = NaN;
        omz(inull) = NaN;
        omr(inull) = NaN;
    end
    
    if size(imA) ~= [nlmesh3,ncmesh3]
        error('Size problem.\n');
    end
    
    if np < 36
        % store images for big panel plot
        oims(i,:,:) = imA;
        mims(i,:,:) = imF;
        rims(i,:,:) = imG;
        cims(i,:,:) = imH;
        tls{i} = titlestr;
    end
    
    % Write phase to binary .pha files and instructions for using them
    %    log_phases(runname, i, imA, imF, imG, imH, imE, uns, mds, urs, ucs);
    %   2012-JUN-25 make everything into meters
%     log_phases(runname, i, imA, imF, imG, imH, imE...
%         , 1.0e3*uns, 1.0e3*mds, 1.0e3*urs, 1.0e3*ucs...
%         , 1.0e3*omr, 1.0e3*omx, 1.0e3*omy, 1.0e3*omz...
%         , 1.0e3*umr, 1.0e3*umx, 1.0e3*umy, 1.0e3*umz);
    write_grids(runname, i, imA, imF, imG, imH, imE...
        , uns, mds, urs, ucs...
        , omr, omx, omy, omz...
        , umr, umx, umy, umz...
        , demx,demy ...
        , OPT.demdescfile, idatatype1);
    
    % set limits of color table
    switch idatatype1 %%TODO handle a different data type for each pair
        case 0
            climits=[-0.5, +0.5];
            if bitget(figopt,1) == 1
                
                ctab = cmapblackzero(1); % black at zero at bottom of color bar
            else
                ctab = colormap('jet');
            end
        otherwise
            climits(1) = nanmin(colvec([imA imB imC imD imE imF imG imH]));
            climits(2) = nanmax(colvec([imA imB imC imD imE imF imG imH]));
            ctab = cmapgraynan;
    end
    
    % set color table
        
%% get centroid of Okada source
%     xcentroid = p1(get_parameter_index('Okada1_Centroid_Easting_in_m____',qnames));
%     ycentroid = p1(get_parameter_index('Okada1_Centroid_Northing_in_m___',qnames));
%     icentroid = NaN;
%     jcentroid = NaN;
%     if isfinite(xcentroid) == 1
%         icentroid = find(abs(yax(isub)-ycentroid) < abs(dy/2.));
%         if numel(icentroid) >= 1
%             icentroid=icentroid(1);
%         end
%     end
%     if isfinite(ycentroid) == 1
%         jcentroid = find(abs(xax(jsub)-xcentroid) < abs(dx/2.));
%         if numel(jcentroid) >= 1
%             jcentroid=jcentroid(1);
%         end
%     end
    
    % draw symbols at centers of sources, UTM coordinates are in meters
%     dotx = xcentroid;
%     doty = ycentroid;  
    dotx = nan;
    doty = nan;
    mysym='w*';
    marksize = 10;
    for ii=1:8
        mysyms{ii} = mysym;
        marksizes(ii) = marksize;
    end
    
    %% set limits of color table
    switch idatatype1 %%TODO handle a different data type for each pair
        case 0
            climits=[-0.5, +0.5];
        otherwise
            % climit(1) = nanmin(nanmin([imA imF imG imH]));
            % climit(2) = nanmax(nanmax([imA imF imG imH]));
            if do_stretch == 1
                climits = quantile(colvec([imA imF imG imH]),[0.05:0.05:0.95]);
            else
                climits(1) = quantile(colvec([imA imF imG imH]),0.05);
                climits(2) = quantile(colvec([imA imF imG imH]),0.95);
            end
    end

    
    %% make 8-panel plot in portrait
    nf=nf+1; h(nf)=utmimage8(imA,imB,imC,imD,imE,imF,imG,imH...
        ,tlA,tlB,tlC,tlD,tlE,tlF,tlG,tlH...
        ,wesn,titlestr,climits,dotx,doty,ctab,1,mysyms,marksizes,idatatype1,datalabel);
    %     ,wesn,titlestr,climit,[Xcorners11 NaN Xcorners21]/1000,[Ycorners11 NaN Ycorners21]/1000,ctab,1,mysyms,marksizes);
    feval(printfun,sprintf('%s_%03d_8PAN',runname,i));
    
    %% make 8-panel plot in landscape
    nf=nf+1; h(nf)=utmimage8landscape(imA,imB,imC,imD,imE,imF,imG,imH...
        ,tlA,tlB,tlC,tlD,tlE,tlF,tlG,tlH...
        ,wesn,titlestr,climits,dotx,doty,ctab,1,mysyms,marksizes,idatatype1,datalabel);
    %     ,wesn,titlestr,climit,[Xcorners11 NaN Xcorners21]/1000,[Ycorners11 NaN Ycorners21]/1000,ctab,1,mysyms,marksizes);
    %feval(printfun,sprintf('%s_%03d_8PANLS.pdf',runname,i),'landscape');
    feval(printfun,sprintf('%s_%03d_8PANLS',runname,i));
    
    %% make 4-panel plot of wrapped phase
    tlA = ''; tlF = '';
    datelabel = '';
    nf=nf+1; h(nf)=utmimage4(imA,imF,imG,imH ...
        ,tlA,tlF,tlG,tlH ...
        ,wesn,titlestr,climits,dotx,doty,ctab,1,mysyms,marksizes...
        ,datelabel,idatatype1,datalabel);
    feval(printfun,sprintf('%s_%03d_4PAN',runname,i));  
     
    
    %% Make profiles
    % decide to show rate or not
    itref = get_parameter_index('Reference_Epoch_in_years________',pnames);
%   20151121 not defined in final estimate
%   if abs(PST.p1(itref)) > 0
%   if abs(PST.p0(itref)) > 0
    if abs(PST1.p0(itref)) > 0
        y0lab='displacement';
        y2lab='mm';
        time_span = 1.0;
    else
        y0lab='displacement rate';
        y2lab='mm/yr';
    end
    y1lab=sprintf('Range %s',y0lab);
    
    
    fprintf(1,'Making E-W profile in range rate for Pair %03d\n',i);
    xt=xmesh3(iprof,:)/1000;
    yt1=uns(iprof,:)*1000/time_span;
    yt2=mds(iprof,:)*1000/time_span;
    %yt2=umr(iprof,:)*1000;
    xlab='Easting (km)';
    %     y1lab='Range change rate (mm/yr)';
    %     y2lab='mm/yr';
    tlab = strcat(sprintf('Northing = %.3f km ',ymesh3(iprof,jprof)/1000),' :',titlestr);
    nf=nf+1;h(nf) = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab);
    %feval(printfun,sprintf('%s_profR_ew_P%03d',runname,i),'landscape');
    feval(printfun,sprintf('%s_profR_ew_P%03d',runname,i));
    
    fprintf(1,'Making N-S profile in range rate for Pair %03d\n',i);
    xt = ymesh3(:,jprof)/1000;
    yt1= uns(:,jprof)*1000/time_span;
    yt2= mds(:,jprof)*1000/time_span;
    xlab='Northing (km)';
    %     y1lab='Range change (mm/yr)';
    %     y2lab='mm/yr';
    
    tlab = strcat(sprintf('Easting = %.3f km ',xmesh3(iprof,jprof)/1000),' :',titlestr);
    nf=nf+1;h(nf) = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab);
    %feval(printfun,sprintf('%s_profR_ns_P%03d',runname,i),'landscape');
    feval(printfun,sprintf('%s_profR_ns_P%03d',runname,i));
    
    if bitget(figopt,3) == 1
        fprintf(1,'Making E-W profiles for Pair %03d\n',i);
        tlab = strcat(sprintf('Northing = %.3f km ',ymesh3(iprof,jprof)/1000),' :',titlestr);
        xt=xmesh3(iprof,:)/1000;
        xlab='Easting (km)';
                
        fprintf(1,'Making E-W profile in east rate for Pair %03d\n',i);
        xt=xmesh3(iprof,:)/1000;
        yt1=omx(iprof,:)*1000/time_span;
        yt2=umx(iprof,:)*1000/time_span;
        xlab='Easting (km)';
        %y1lab='Eastward velocity (mm/yr)';
        y1lab=sprintf('Eastward %s',y0lab);
        
        %         y2lab='mm/yr';
        nf=nf+1;h(nf) = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab);
        feval(printfun,sprintf('%s_profE_ew_P%03d',runname,i));
        
        fprintf(1,'Making E-W profile in north for Pair %03d\n',i);
        yt1=omy(iprof,:)*1000/time_span;
        yt2=umy(iprof,:)*1000/time_span;
        %y1lab='Northward velocity (mm/yr)';
        y1lab=sprintf('Northward %s',y0lab);
        
        %         y2lab='mm/yr';
        nf=nf+1;h(nf) = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab);
        feval(printfun,sprintf('%s_profN_ew_P%03d',runname,i));
        
        fprintf(1,'Making E-W profile in vertical rate for Pair %03d\n',i);
        yt1=omz(iprof,:)*1000/time_span;
        yt2=umz(iprof,:)*1000/time_span;
        %y1lab='Upward velocity (mm/yr)';
        y1lab=sprintf('Upward %s',y0lab);
        
        %         y2lab='mm/yr';
        nf=nf+1;h(nf) = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab);
        feval(printfun,sprintf('%s_profU_ew_P%03d',runname,i));
        
        %   2012-10-04 values in meters
        kq = kq+1;
        q0(kq)     = NaN;
        q1(kq)     = nanmin(colvec(umz));
        qsig(kq)   = nanstd(colvec(umz));
        qnames{kq} = sprintf('Pair_%05d_MiniVerti_m_MOD__',i);
        kq = kq+1;
        q0(kq)     = NaN;
        q1(kq)     = nanmax(colvec(umz));
        qsig(kq)   = nanstd(colvec(umz));
        qnames{kq} = sprintf('Pair_%05d_MaxiVerti_m_MOD__',i);
        
        %         %   2012-10-22 values in meters
        %         iq = iq+1;
        %         q0(iq)     = NaN;
        %         q1(iq)     = umx(iprof,jprof);
        %         qsig(iq)   = NaN;
        %         qnames{iq} = sprintf('Pair_%05d_DisXatCen_m_MOD__',i);
        %         iq = iq+1;
        %         q0(iq)     = NaN;
        %         q1(iq)     = umy(iprof,jprof);
        %         qsig(iq)   = NaN;
        %         qnames{iq} = sprintf('Pair_%05d_DisYatCen_m_MOD__',i);
        %         iq = iq+1;
        %         q0(iq)     = NaN;
        %         q1(iq)     = umz(iprof,jprof);
        %         qsig(iq)   = NaN;
        %         qnames{iq} = sprintf('Pair_%05d_DisZatCen_m_MOD__',i);
        
        %   2012-10-25 values in meters
%         if numel(icentroid) == 1 && numel(jcentroid) == 1
%             if isfinite(icentroid) == 1 && isfinite(jcentroid) == 1
%                 kq = kq+1;
%                 q0(kq)     = NaN;
%                 q1(kq)     = umx(icentroid,jcentroid);
%                 qsig(kq)   = NaN;
%                 qnames{kq} = sprintf('Pair_%05d_ModDisXCentroid_m',i);
%                 kq = kq+1;
%                 q0(kq)     = NaN;
%                 q1(kq)     = umy(icentroid,jcentroid);
%                 qsig(kq)   = NaN;
%                 qnames{kq} = sprintf('Pair_%05d_ModDisYCentroid_m',i);
%                 kq = kq+1;
%                 q0(kq)     = NaN;
%                 q1(kq)     = umz(icentroid,jcentroid);
%                 qsig(kq)   = NaN;
%                 qnames{kq} = sprintf('Pair_%05d_ModDisZCentroid_m',i);
%             end
%         end
    end
    
    % calculate some derived parameters
    %     iq=iq+1;
    %     q0(iq) = max(max(abs(imA)));
    %     q1(iq) = max(max(abs(imE)));
    %     qsig(iq) = NaN;
    %     qnames{iq} = sprintf('Pair_%05d_MaxAbs_Cycles_OBS',i);
    %iq=iq+1;
    %     q0(iq) = min(min(mdl1));
    %     q1(iq) = max(max(mdl1));
    %   2012-06-13 values in millimeters
    %     q0(iq) = 1.0e3*DST.mpercy(1)*double(min(min(mdl1)))/DNPC;
    %     q1(iq) = 1.0e3*DST.mpercy(1)*double(max(max(mdl1)))/DNPC;
    %     %   2012-06-26 values in meters
    %     q0(iq) = DST.mpercy(1)*double(min(min(mdl1)))/DNPC;
    %     q1(iq) = DST.mpercy(1)*double(max(max(mdl1)))/DNPC;
    %     qsig(iq) = NaN;
    %     qnames{iq} = sprintf('Pair_%05d_MinMaxRng_m_MOD__',i);
    %     iq=iq+1;
    %     q0(iq) = min(min(urs));
    %     q1(iq) = max(max(urs));
    %     qsig(iq) = NaN;
    %     qnames{iq} = sprintf('Pair_%05d_MinMaxRng_m_RES__',i);
    %     iq=iq+1;
    %     q0(iq) = min(min(ucs));
    %     q1(iq) = max(max(ucs));
    %     qsig(iq) = NaN;
    %     qnames{iq} = sprintf('Pair_%05d_MinMaxRng_m_DEV__',i);
    %     iq=iq+1;
    %     q0(iq) = uns(iprof,jprof);
    %     q1(iq) = mds(iprof,jprof);
    %     qsig(iq) = NaN;
    %     qnames{iq} = sprintf('Pair_%05d_Unwrap_Obs_Mod_mm',i);
    %     mds(iprof,jprof)/1.0e3;
    %   2012-10-04 values in meters
%     iq = iq+1;
%     q0(iq)     = DST.mpercy(1)*nanmin(colvec(double(mdl0)))/DNPC;
%     q1(iq)     = DST.mpercy(1)*nanmin(colvec(double(mdl1)))/DNPC;
%     qsig(iq)   = DST.mpercy(1)*nanstd(colvec(double(mdl1)))/DNPC;
%     qnames{iq} = sprintf('Pair_%05d_MiniRange_m_MOD__',i);
%     iq = iq+1;
%     q0(iq)     = DST.mpercy(1)*nanmax(colvec(double(mdl0)))/DNPC;
%     q1(iq)     = DST.mpercy(1)*nanmax(colvec(double(mdl1)))/DNPC;
%     qsig(iq)   = DST.mpercy(1)*nanstd(colvec(double(mdl1)))/DNPC;
%     qnames{iq} = sprintf('Pair_%05d_MaxiRange_m_MOD__',i);
%     
%     
%     %  PRINT OUT THE DERIVED PARAMETERS
%     qnames = truncate_parameter_names(qnames);
%     iq1 = iq2+1;
%     iq2 = iq;
%     uqb(iq1:iq2)=NaN;
%     lqb(iq1:iq2)=NaN;
%     
%     for j=iq1:iq2
%         qflags{j} = 'D#';
%         adj = q1(j)-q0(j);
%         
%         outfmt = getfmt(q1(j),qnames{j});
%         
%         fprintf(1        ,outfmt,qflags{j},j+mparam,qnames{j} ,q0(j),q1(j),adj,qsig(j),sadj,(uqb(j)-lqb(j))/2.0);
%         fprintf(fidtxtout,outfmt,qflags{j},j+mparam,qnames{j}, q0(j),q1(j),adj,qsig(j),sadj,(uqb(j)-lqb(j))/2.0);
%     end
    
%     % print corners of Okada models
%     iii=get_parameter_index('Okada1_Length_in_m______________',pnames);
%     if p0(iii) > 0
%         fprintf(1,'Initial UTM coordinates of 4 corners, upper center, and Centroid of Okada1\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},Xcorners10(iiii),Ycorners10(iiii),Hcorners10(iiii));
%         end
%         fprintf(1,'Final   UTM coordinates of 4 corners, upper center, and Centroid of Okada1\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},Xcorners11(iiii),Ycorners11(iiii),Hcorners11(iiii));
%         end
%         fprintf(1,'Initial Lon Lat of 4 corners, upper center, Centroid of Okada1\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},LonCorners10(iiii),LatCorners10(iiii),Hcorners10(iiii));
%         end
%         fprintf(1,'Final   Lon Lat of 4 corners, upper center, Centroid of Okada1\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},LonCorners11(iiii),LatCorners11(iiii),Hcorners11(iiii));
%         end
%     end
%     iii=get_parameter_index('Okada2_Length_in_m______________',pnames);
%     if p0(iii) > 0
%         fprintf(1,'Initial UTM coordinates of 4 corners, upper center, and Centroid of Okada2\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},Xcorners20(iiii),Ycorners20(iiii),Hcorners20(iiii));
%         end
%         fprintf(1,'Final   UTM coordinates of 4 corners, upper center, and Centroid of Okada2\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},Xcorners21(iiii),Ycorners21(iiii),Hcorners21(iiii));
%         end
%         fprintf(1,'Initial Lon Lat of 4 corners, upper center, Centroid of Okada2\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},LonCorners20(iiii),LatCorners20(iiii),Hcorners20(iiii));
%         end
%         fprintf(1,'Final   Lon Lat of 4 corners, upper center, Centroid of Okada2\n');
%         for iiii=[1 2 3 4 7 10]
%             fprintf(1,'%s %12.4f %12.4f %12.4f\n',Ncorners{iiii},LonCorners21(iiii),LatCorners21(iiii),Hcorners21(iiii));
%         end
%     end
    
    
    % do not open too many windows
    if np > 2
        close all
    end
    
end % loop over pairs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make multi-panel plots
% commented out 11/14/11 to allow ensemble to run
if np > 1 && np < 36
    %wesn = [lonmin lonmax latmin latmax];
    wesn = [xmin xmax ymin ymax];
    
    % pad out the rest of the panels
    for i=np+1:nfull
        oims(i,:,:) = zeros(nlmesh3,ncmesh3);
        mims(i,:,:) = zeros(nlmesh3,ncmesh3);
        rims(i,:,:) = zeros(nlmesh3,ncmesh3);
        cims(i,:,:) = zeros(nlmesh3,ncmesh3);
        tls{i} = 'Models';
    end
    % utmmageNM(ims,tls,nrows,mcols....
    %             ,wesn,titlestr,climit,dotutmx,dotutmy,ctab)
    
    climits=[-0.5, +0.5];
    
    %nf=nf+1;h(nf)=figure;
    utmimageNM(oims,tls,panelrows,panelcols,wesn...
        ,'Observed Phase Values',climits,0.,0.,ctab)
    feval(printfun,sprintf('%s_ALLOBS',runname));
    
    %nf=nf+1;h(nf)=figure;
    utmimageNM(mims,tls,panelrows,panelcols,wesn...
        ,'Final Modeled Phase Values',climits,0,0,ctab)
    feval(printfun,sprintf('%s_ALLMOD',runname));
    
    %nf=nf+1;h(nf)=figure;
    utmimageNM(rims,tls,panelrows,panelcols,wesn...
        ,'Final Residual Phase Values',climits,0,0,ctab)
    feval(printfun,sprintf('%s_ALLRES',runname));
    
    %nf=nf+1;h(nf)=figure;
    utmimageNM(cims,tls,panelrows,panelcols,wesn...
        ,'Final Deviations in Phase',climits,0,0,ctab)
    feval(printfun,sprintf('%s_ALLDEV',runname));
    
end


%save ; % do not save all these variables becauase they overwrite some others

fclose(fidtxtout);

% % append derived parameters (labeled q) to estimated parameters
% pqnum = mparam + iq;
% pq0 = [colvec(p0); colvec(q0)];
% pq1 = [colvec(p1); colvec(q1)];
% %qsig = nan(iq,1);
% pqsig = [colvec(psig); colvec(qsig)];
% pqnames = [pnames qnames];
% qbounds = nan(iq,2);
% pqbounds = [bounds; qbounds];
% pqflags = [pflags qflags];
% qscl = nan(iq,1);
% pqscl = [colvec(pscl); colvec(qscl)];

% write the parameter structure to a file
%PSTF = build_pst(fitfun,pqnum,pq0,pq1,pqsig,pqnames,pqbounds,datafilename,pqscl,pqflags,timefun);
% write_pst(PSTF,'PST.OUT');
% ierr = print_parameters_nicely(PSTF,fnparout);
clear h;

%% This file is too big to save.
%% save('gipht.mat');
fdelete('gipht.mat');

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));

return



