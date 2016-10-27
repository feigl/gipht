% gipht_step1
% read input files
% define subregion
% read phase data
% select pixels
% 2009-JUL-14 Big clean up
% 2009-DEC-08 Allow pselect = 7
% 2010-NOV-04 Write 2nd DEM descriptor for subregion
% 2012-JAN-12 try to pass test on simulated data
% 2012-OCT-01 back to filtered psp instead of unfiltered pha
% 20160524    read GMT grid files from GMT5SAR

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

% [xcenter, ycenter, halfwidth, halfheight, npix, pselect, tquake...
%     , unitv0, ithresh, maxcmd, pixinpatch, maxpix, ianneal, nprocessors, interpcell ...
%     , ilist, fnparin, fnparout, objfun, fitfun, demdescfile, orbfile, cohfile...
%     , mpercy, datafilename, nsaruns, parmfilename, saopt6, figopt, printfun, orbopt...
%     , pha2qlsname, phaseprefix, surrogate, verbose]...
%     = read_input_controls('gipht.in',runname);

if fexist('gipht.mat') == 1
    load('gipht.mat');
end

if exist('OPT','var') == 0
    OPT = struct([]);
end

%% unpack options
OPT = read_input_controls(OPT);
% functions
%objfun       = OPT.objfun;
fitfun       = OPT.fitfun;
fitfun_exact = OPT.fitfun;
printfun     = OPT.printfun;
timefun      = OPT.timefun;
% variables
figopt       = OPT.figopt;
ianneal      = OPT.ianneal;
ilist        = OPT.ilist;
pselect      = OPT.pselect;
npix         = OPT.npix;
orbopt       = OPT.orbopt;
unitv0       = OPT.unitv;
saopt6       = OPT.saopt6;
% file names
datafilename = OPT.datafilename;
fnparin      = OPT.fnparin;
fnparout     = OPT.fnparout;
fntxtout     = OPT.fnsumout;
demdescfile  = OPT.demdescfile;
orbfile      = OPT.orbfile;
nsaruns      = OPT.nsaruns;
verbose      = OPT.verbose;
nprocessors  = OPT.nprocessors;


% how to handle statistics
istatcode = 0;

% name for DEM descriptor file that describes subregion
%demdescfile2 = sprintf('gipht_subregion.dat');

% name for DST file
%fnamedst = 'dst_sample.dst';

% name for PST file
%fnamepstin = 'pst_in.pst';

% name of file to save orbits
orbsavefile = sprintf('orbits.mat');
orbits_loaded = 0;

% open some output files
fidtxtout=fopen(OPT.fnsumout,'w');
if fidtxtout <= 0
    error(sprintf('Cannot open output file %s\n',OPT.fnparout));
end
fprintf(fidtxtout,'%s %s %s\n',OPT.fnparout,runname,versionstr);
%fprintf(fidtxtout,'%80s\n',splashtext);

% % deal with parallel processing if requested
% if nprocessors > 1
%     if numel(which('matlabpool')) > 0
%         %         fprintf(1,'Closing any lingering workers in matlab pool.\n');
%         %         matlabpool close force
%         fprintf(1,'Attempting to open matlabpool with %d processors\n',nprocessors);
%         matlabpool('open',nprocessors)
%         if matlabpool('size') == nprocessors
%             fprintf(1,'Opened matlabpool with %d processors\n',nprocessors);
%         else
%             warning('Could not open matlab pool for distributed computing');
%             nprocessors = 1;
%         end
%     else
%         warning(sprintf('Request is for nprocessors = %d BUT matlab distributed tool kit not installed.\n',nprocessors));
%         nprocessors = 1;
%     end
% end

% % directory to look for orbit files
% if numel(orbfile) > 0
%     if orbfile(end) == '/'
%         orbdir = orbfile;
%     else
%         orbdir = '';
%     end
% else
%     orbdir = '';
% end

%% list of file names
if fexist(OPT.ilist) > 0
    % read the phase file names from file_names.dat
    disp('Name of list of phase files');
    [pfnames, mdate, imast, sdate, islav, hamb, ddays, t1, t2, idatatypes, mpercys] = read_file_names(OPT.ilist);
else
    error(sprintf('Could not find file named %s listing phase files\n',OPT.ilist));
end

% Check that all files exist
nk=0;
for i=1:numel(pfnames)
    if fexist(pfnames{i}) > 0
        nk=nk+1;
    else
        warning(sprintf('Cannot find  %80s\n',pfnames{i}));
    end
end
nk


%% select data type and objective function
ndatatypes = sum(idatatypes)/nk;
switch ndatatypes
    case 0  % wrapped phase
        idatatype1 = 0;
        FACTIN = 1; % .grd file contains radians
        DNPC = 2 * pi;     
        datalabel = '[cycles]';
        objfun = 'funcostrarc';        % Objective Function mininum angle,  assumes zero mean, using arc function in radians
        objlabel = '[cycles]';
    case -1; % east component of gradient
        idatatype1 = -1;
        FACTIN = 1; % grd file contains dimensionless strain
        DNPC = 1;
        datalabel = '[dimless]';
        objfun = 'funcoststdnres';      % Objective function is sample standard deviation of normalized residual (should equal sqrt(chi2))
        objlabel = '[dimless]';
    case 2   % range change in meters after unwrapping
        idatatype1 = 2;
        FACTIN = 1.; % grd file contains meters
        DNPC = 1.0e-3;     
        datalabel = '[mm]';
        objfun = 'funcoststdnres';      % Objective function is sample standard deviation of normalized residual (should equal sqrt(chi2))
        objlabel = '[dimless]';
    otherwise
        error(sprintf('unknown ndatatypes %d\n',ndatatypes));
end


%% count pairs
if nk == numel(pfnames)
    %np = length(imast)
    np = nk;
    fprintf(1,'Number of interferogram files (pairs) to adjust = %d\n',np);
else
    error(sprintf('Missing %d of %d interferogram files\n'...
        ,numel(pfnames)-nk,numel(pfnames)));
end


% %% handle orbits
% if numel(orbfile) > 0
%     for i=1:numel(imast);
%         % Accomodate numbering scheme for Okmok used by Zhong Lu
%         jmast =  mod(abs(imast(i)),100000);
%         jslav =  mod(abs(islav(i)),100000);
%         %pfnames{i} = sprintf('../In%4d_%4d/psp_%4d_%4d_ort.pha',jmast,jslav,jmast,jslav);
%         temponame=sprintf('%s%d.orb',orbdir,jmast);
%         if fexist(temponame) > 0
%             ofnames1{i} = temponame;
%         else
%             fprintf(1,'WARNING: Could NOT find orbit file named %s. Substituting %s\n',temponame,orbfile);
%             ofnames1{i} = orbfile;
%         end
%         temponame=sprintf('%s%d.orb',orbdir,jslav);
%         if fexist(temponame) > 0
%             ofnames2{i} = temponame;
%         else
%             fprintf(1,'WARNING: Could NOT find orbit file named %s. Substituting %s\n',temponame,orbfile);
%             ofnames2{i} = orbfile;
%         end
%     end
%     for i=1:numel(ofnames1)
%         if fexist(ofnames1{i}) > 0
%             fprintf(1,'Master orbit file for pair %d is named %s\n',i,ofnames1{i});
%         else
%             error(sprintf('Cannot find  %80s\n',ofnames1{i}));
%         end
%     end
%     for i=1:numel(ofnames2)
%         if fexist(ofnames2{i}) > 0
%             fprintf(1,'Slave orbit file for pair %d is named %s\n',i,ofnames2{i});
%         else
%             error(sprintf('Cannot find  %80s\n',ofnames2{i}));
%         end
%     end
% end


%% approximate baseline in meters from Ha in meters
bperp = ones(size(imast));

%% find the disconnected trees
[trees, DD, tepochs, iepochs, iuniqorbs, uniqdates] = findtrees2(t1,t2);

%% number of epochs
me = length(tepochs);
fprintf(1,'Number of distinct epochs overall me = %d\n',me);
mfam = me - rank(DD);

%% number of trees
[ntrees,ms] = size(trees);
if ntrees == mfam
    fprintf(1,'Number of distinct trees (components) overall: ntrees = %d\n',ntrees);
else
    error(sprintf('Error counting components: %d %d\n',ntrees,mfam));
end

%% loop over trees to count the number of pairs needed
np4 = 0;
for i=1:ntrees
    fprintf(1,        'Species %1s:',char(i+64));
    fprintf(fidtxtout,'Species %1s:',char(i+64));
    for j=1:ms
        if isfinite(trees(i,j)) == 1
            fprintf(1,         ' %5d',trees(i,j));
            fprintf(fidtxtout, ' %5d',trees(i,j));
            if j > 1
                np4 = np4 +1;
            end
        end
    end
    fprintf(1,        '\n');
    fprintf(fidtxtout,'\n');
end

%% loop over epochs to list pairs
np2 = 0;
np3=0;
for k=1:me
    % index of master
    %kmast = find(DD(k,:) == -1);
    % index of pair
    kp = find(DD(:,k) == -1); % epoch is master
    if numel(kp) == 0
        kp = find(DD(:,k) == +1); % epoch is slave
    end
    if numel(kp) > 0
        np2 = np2 + 1;
    else
        error('pair index not found!');
    end
    
    [is,jm] = find(trees == k);
    % is is index to trees
    % jm is index of member within trees
    for kk=1:numel(kp)
        np3 = np3+1;
        tmp=sprintf('Species %1s Member %2d Pair %3d M %6d S %6d %80s',char(is+64),jm,kp(kk),imast(kp(kk)),islav(kp(kk)),char(pfnames{kp(kk)}));
        pairlist{np3}=tmp;
        %fprintf(1,'%6d %6d %80s\n',imast(kp(kk)),islav(kp(kk)),pfnames{kp(kk)});
    end
end
pairlist = sort(pairlist);
for k=1:numel(pairlist)
    fprintf(1,        '%s\n',char(pairlist{k}));
    fprintf(fidtxtout,'%s\n',char(pairlist{k}));
end

%% TODO replace with Elena's functions
% estimate the pseudoabsolute baselines
bpest = adjustbp(tepochs,DD,bperp, trees,iuniqorbs, uniqdates);
nf=nf+1;h(nf)=figure;
%plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, 0, 'orbital separation (Bperp) [m]');
%ktour = plotbp(tepochs, bpest, DD, trees,iuniqorbs, uniqdates, 3);

% Call Elena's function 20150901
% function h = plot_trees(tepochs, scores, DD, trees,xlab,ylab)
plot_trees(tepochs, bpest, DD, trees, 'year', 'B_perp [m]');
feval(printfun,sprintf('%s_TREES',runname));

if np > np4
    warning(sprintf('Number of pairs (%d) is larger than necessary (%d). Check for redundant pairs in list above.\n',np,np4));
end

%% should estimate the pseuodabsolute doppler values
dops = zeros(size(tepochs));

%% find out about the DEM
[isgeo,y1,x1,nl,nc,l1,c1,ml,mc,dl,dc,fi2,lat0,lon0,y0,x0,hemisphere,iutmzone,isgmtgrid] = read_dem_descriptor(demdescfile);

disp('number of columns in each interferogram');ncols = mc
disp('number of lines in each interferogram '); nrows = ml

[demx,demy,demz] = grdread3(demdescfile);
nf=nf+1;h(nf) = map_grd(demdescfile,'copper');
xsubmin = nanmin(demx);
xsubmax = nanmax(demx);
ysubmin = nanmin(demy);
ysubmax = nanmax(demy);
xcenter = mean([xsubmin, xsubmax]);
ycenter = mean([ysubmin, ysubmax]);


%% Decide how to select pixels
switch pselect
    case 0
        fprintf(1,'Selecting all pixels.\n');
    case 1
        fprintf(1,'Selecting pixels randomly.\n');
    case 2
        fprintf(1,'Reading previous pixel indices from existing files ikeep.mat and jkeep.mat.\n');
    case 3
        fprintf(1,'Selecting pixels using quadtree subsampling, then randomly\n');
%     case 4
%         fprintf(1,'Selecting pixels based on coherence.\n');
    case 5
        fprintf(1,'Selecting pixels using pha2qls resampling.\n');
%     case 6
%         fprintf(1,'Selecting pixels based on quadtree of stacked phase.\n');
    case 7
        fprintf(1,'Selecting pixels using pha2qls. Values are range gradient.\n');
    case 9
        fprintf(1,'Selecting pixels randomly. Values are range gradient.\n');
    otherwise,
        error(sprintf('Unkown value of pselect (%d)\n',pselect));
end


%% set pointers for indices of pixels
ippix1 = zeros(np,1); i1 = 0; % index to first pixel in pair
ippix2 = zeros(np,1); i2 = 0; % index to last pixel in pair

%% loop over pairs
nk=0;      % count the number of pairs 
ndata1=0;  % count the number of data in the current pair
i1=1;      % index to start of data vector
i2=0;      % index to end of data vector
%kk=0;      % count the number of data
kindex = zeros(np,1);
kmasts = zeros(np,1);
kslavs = zeros(np,1);

for i = 1:np
    %% index master and slave
    kmast = find(DD(i,:) == -1);
    kslav = find(DD(i,:) == +1);
    
    %%read phase data from interferogram
    %2009-JUN-18 phaimg = read_pha(fn0,ncols)/256; % in cycles
    
    fn0 = pfnames{i};
    nbytes = fsize(fn0);
    if  nbytes > 0
        fprintf(1,'Extracting information from GMT UTM grid file named %s\n',fn0);
        INFO = grdinfo3(fn0)
        nf=nf+1;h(nf) = map_grd(fn0,cmapblackzero(1));
        feval(printfun,sprintf('%s_obs_P%02d',runname,i));  
        [grdx,grdy,grdd] = grdread3(fn0);       
     else
        phaimg = [];
        warning(sprintf('Phase file named %s is non-existant or empty\n',fn0));
        nerrors = nerrors + 1;
    end
    
    %% decide about masking
    switch pselect
        case {5,7}
            fprintf(1,'Masking.\n');
            fnmsk = strrep(fn0,'qphase','qmaskn');
            if fexist(fnmsk) == 1
                fprintf(1,'Opening mask file named %s ...\n',fnmsk);
                [mskx,msky,mskd] = grdread3(fnmsk);
            else
                warning(sprintf('Could not open mask file named %s . Taking all pixels',fnmsk));
                mskx = grdx;
                msky = grdy;
                mskd = grdd;
            end
        otherwise
            warning(sprintf('Not masking.\n'));
            mskx = grdx;
            msky = grdy;
            mskd = grdd;
    end

    
    %% verify that grids are same size
    if numel(grdx) ~= ncols || numel(mskx) ~= ncols
        error(sprintf('number of columns in data file (%d) or mask file (%d) is not equal to the number of columns (%d) in DEM\n'...
            ,numel(grdx),numel(mskx),ncols));  
    end  
    if numel(grdy) ~= nrows || numel(msky) ~= nrows
        error(sprintf('number of rows in data file (%d) or mask file (%d) is not equal to the number of columns (%d) in DEM\n'...
            ,numel(grdy),numel(msky),nrows));  
    end  
    
    %% select the pixels
    kok = select_pixels(pselect,npix,grdd,mskd);
    [iok,jok] = ind2sub(size(grdd),kok);
    
    %% extract the data and load into arrays
    ndata1 = numel(kok);
    i2 = i2+ndata1;
    phao(i1:i2)   = double(grdd(kok))*FACTIN;
    %% TODO fix the next two lines X and Y coordinates from vectors
    xyzm(1,i1:i2) = double(grdx(jok));
    xyzm(2,i1:i2) = double(grdy(iok));
    % elevation is a grid
    xyzm(3,i1:i2) = double(demz(kok));
    alat(i1:i2) = nan(size(kok));
    alon(i1:i2) = nan(size(kok));
    kindex(i1:i2) = i;
    kmasts(i1:i2) = kmast*ones(size(kok));
    kslavs(i1:i2) = kslav*ones(size(kok));
    
    fprintf(1,'Minimal values in meters         for xyzm %10.1f %10.1f %10.1f \n',min(xyzm(1,i1:i2)),min(xyzm(2,i1:i2)),min(xyzm(3,i1:i2)));
    fprintf(1,'Maximal values in meters         for xyzm %10.1f %10.1f %10.1f \n',max(xyzm(1,i1:i2)),max(xyzm(2,i1:i2)),max(xyzm(3,i1:i2)));
    
    %% make a title for subsequent plots
    titlestr = sprintf('Phase: Pair %3d %s orbs %5d %5d NROWS = %d by NCOLS= %d'...
        ,i,strrep(fn0,'_','\_'),iuniqorbs(kmast),iuniqorbs(kslav),nrows,ncols);

    %% plot wrapped phase in histogram
    nf=nf+1;h(nf)=figure;title(titlestr);
    hist2(phao(i1:i2)/DNPC,64);hold on;
    axis xy;xlabel(datalabel); ylabel('number of occurences');
    feval(printfun,sprintf('%s_Histogram_P%02d',runname,i));
    
    %% plot data as a function of topographic elevation
    nf=nf+1;h(nf)=figure;hold on;title(titlestr);
    subplot(4,1,1);plot(xyzm(1,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k.');xlabel('Easting (km)');ylabel(datalabel);
    subplot(4,1,2);plot(xyzm(2,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k.');xlabel('Northing (km)');ylabel(datalabel);
    subplot(4,1,3);plot(xyzm(3,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k.');xlabel('topographic elevation (km)');ylabel(datalabel);
    subplot(4,1,4);rdist = sqrt((xyzm(1,i1:i2)-xcenter).^2 +(xyzm(2,i1:i2)-ycenter).^2);
    plot(rdist/1e3, double(phao(i1:i2))/DNPC,'k.'); xlabel('Radial distance from center of subregion (km)');ylabel(datalabel);
    feval(printfun,sprintf('%s_PositionVsPhase_P%02d',runname,i));
   
    %% deal with orbits
    [unitv, orbvm, orbvs] = orbit_handler(i,i1,i2,orbfile, orbits_loaded, orbopt, unitv0);
 
    % close some windows if too many
    if np > 2
        close all
    end
    
   %% update index pointing to end of data vectors
   i1 = i1+ndata1;  
   %% update count of pairs
   nk = nk+1;
end % loop over files

% % Check that orbits are consistent with data
% [ndum,nunitv]=size(unitv);
% if ndum == 3 && nunitv == ndata1
%     orbits_loaded = 1;
%     fprintf(1,'Saving orbit information for future use in file named %s\n',orbsavefile);
%     if exist('orbvm','var') ~= 1
%         orbvm = zeros(6,ndata1);
%     end
%     if exist('orbvs','var') ~= 1
%         orbvs = zeros(6,ndata1);
%     end
%     save(orbsavefile,'alat','alon','unitv','orbvm','orbvs');
% else
%     orbits_loaded = 0;
%     fprintf(1,'Dimension mismatch number of unit vectors (%d) different from number of pixels(%d)\n',nunitv,ndata1);
%     fprintf(1,'Deleting save file %s\n',orbsavefile);
%     fdelete(orbsavefile);
%     fprintf(1,'Please run again.\n');
%     error(sprintf('Dimension mismatch number of unit vectors (%d) different from number of pixels(%d)\n',nunitv,ndata1));
% end
% 


disp('Number of pairs with data');np = nk

disp('Number of data');ndata = numel(phao)
disp('Number of data');ndata1

if ndata ~= ndata1
    error(sprintf('NDATA = %d but NDATA1 = %d\n',ndata,ndata1));
end


% i1 = 1;
% i2 = ndata;

%% build DST structure with data
% 2011-OCT-03 NOTE!!! phao is written in RADIANS to DST file

%% Step sizes equal differences in coordinates
dx = grdx(2)-grdx(1);
dy = grdy(2)-grdy(1);
dz(1) = 0;
dz(2:ndata) = colvec(diff(xyzm(3,:)));

% do not need quad-tree indices
%if ismember(pselect,[5,7]) == 0
    qii1 = zeros(ndata,1);
    qii2 = zeros(ndata,1);
    qjj1 = zeros(ndata,1);
    qjj2 = zeros(ndata,1);
%end

% measurement uncertainty
phasig = ones(size(phao));
DST=build_dst(fitfun,xyzm,tepochs,bpest,dops,DD,unitv...
    ,phao,NaN*ones(size(phao)),ippix1,mpercys,idatatype1...
    ,dx,dy,dz,orbvm,orbvs,alon,alat...
    ,qii1,qii2,qjj1,qjj2...
    ,phasig,kindex,kmasts,kslavs);
% no need to write structure
%ierr = write_dst(DST,fitfun,fnamedst);
%check_struct(DST,0); % do not print
check_struct(DST,1); % print min, max, values


% % % plot the pseudoabsolute baselines
% % nf=nf+1;h(nf)=figure;
% % plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, 0);
% % ktour = plotbp(tepochs, bpest, DD, trees,iuniqorbs, uniqdates, 3);
% % feval(printfun,sprintf('%s_%02d',mfilename,nf));
%

% %
% figure;plot(alon,unitv(1,:),'r+');ylabel('east component of unit vector');axis tight
% figure;plot(alon,unitv(2,:),'r+');ylabel('north component of unit vector');axis tight
% figure;plot(alon,unitv(3,:),'r+');ylabel('up component of unit vector');axis tight
%
% save unitv unitv

clear inan imout iok i2dem isort phaimg phao rng xyzm;
clear qgx qgy qi1 qi2 qim qim2 qj1 qj2 qjm qjm2 qqp qqq qsp qkw qlist;
clear qi12 qiim qj12 qjjm;
clear h;
save('gipht.mat');

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));
