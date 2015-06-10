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

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

[xcenter, ycenter, halfwidth, halfheight, npix, pselect, tquake...
    , unitv0, ithresh, maxcmd, pixinpatch, maxpix, ianneal, nprocessors, interpcell ...
    , ilist, txtinname, txtoutname, objfun, fitfun, demdescfile, orbfile, cohfile...
    , mpercy, datafilename, nsaruns, parmfilename, saopt6, figopt, printfun, orbopt...
    , pha2qlsname, phaseprefix, surrogate, verbose]...
    = read_input_controls('gipht.in',runname);

% record fitting function
fitfun_exact = fitfun; 

% how to handle statistics
istatcode = 0;

% name for DEM descriptor file that describes subregion
demdescfile2 = sprintf('gipht_subregion.dat');

% name for DST file
fnamedst = 'dst_sample.dst';

% name for PST file
fnamepstin = 'pst_in.pst';

% name of file to save orbits
orbsavefile = sprintf('orbits.mat');
orbits_loaded = 0;

% choose data type
if ismember(pselect,[7,9])
    idatatype1 = -1; % gradients
else
    idatatype1 = 0;  % wrapped phase
end

disp('Scale factor for phase in DN per cycle:')
% if numel(findstr(objfun,'16')) > 0
%     DNPC = 256^2;  % for gradients
% end
DNPC = 2 * pi


% open some output files
fidtxtout=fopen(txtoutname,'w');
if fidtxtout <= 0
    error(sprintf('Cannot open output file %s\n',txtoutname));
end
fprintf(fidtxtout,'%s %s %s\n',txtoutname,runname,versionstr);
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

% directory to look for orbit files
if numel(orbfile) > 0
    if orbfile(end) == '/'
        orbdir = orbfile;
    else
        orbdir = '';
    end
else
    orbdir = '';
end

% list of file names
if strfind(ilist,'file_names') > 0
    % read the phase file names from file_names.dat
    disp('Name of list of phase files');
    [pfnames, mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_file_names(ilist);
else
    % read the list of interferograms to get pairs marked with little 'a'
    disp('Name of list of interferograms');
    [mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_igram_list(ilist,'a');
    
    % define the file names based on usual convention using DIAPASON and DTOOLS
    for i=1:numel(imast);
        jmast =  abs(imast(i));
        jslav =  abs(islav(i));
        dirname = sprintf('../IN');
        if dexist(dirname) == 0
            dirname = sprintf('../In%4d_%4d',jmast,jslav);
        end
        %       fn0 = sprintf('psp_%4d_%4d_ort.pha',jmast,jslav);
        %       20120926 use unfiltered phase values
        %       fn0 = sprintf('pha_%4d_%4d_ort.pha',jmast,jslav);
        %       2012-OCT-01 back to filtered psp instead of unfiltered pha
        %       fn0 = sprintf('psp_%4d_%4d_ort.pha',jmast,jslav);
        %       2012-OCT-04 use prefix
        fn0 = sprintf('%3s_%4d_%4d_ort.pha',phaseprefix,jmast,jslav);
        if fexist(fn0) == 1
            pfnames{i} = fn0;
        else
            if dexist(dirname) == 1
                pfnames{i} = sprintf('%s/%s',dirname,fn0);
            else
                warning(sprintf('Cannot find directory named %s\n',dirname));
                pfnames{i} = sprintf('%s/%s','./',fn0);
            end
        end
    end
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

if numel(orbfile) > 0
    for i=1:numel(imast);
        % Accomodate numbering scheme for Okmok used by Zhong Lu
        jmast =  mod(abs(imast(i)),100000);
        jslav =  mod(abs(islav(i)),100000);
        %pfnames{i} = sprintf('../In%4d_%4d/psp_%4d_%4d_ort.pha',jmast,jslav,jmast,jslav);
        temponame=sprintf('%s%d.orb',orbdir,jmast);
        if fexist(temponame) > 0
            ofnames1{i} = temponame;
        else
            fprintf(1,'WARNING: Could NOT find orbit file named %s. Substituting %s\n',temponame,orbfile);
            ofnames1{i} = orbfile;
        end
        temponame=sprintf('%s%d.orb',orbdir,jslav);
        if fexist(temponame) > 0
            ofnames2{i} = temponame;
        else
            fprintf(1,'WARNING: Could NOT find orbit file named %s. Substituting %s\n',temponame,orbfile);
            ofnames2{i} = orbfile;
        end
    end
    for i=1:numel(ofnames1)
        if fexist(ofnames1{i}) > 0
            fprintf(1,'Master orbit file for pair %d is named %s\n',i,ofnames1{i});
        else
            error(sprintf('Cannot find  %80s\n',ofnames1{i}));
        end
    end
    for i=1:numel(ofnames2)
        if fexist(ofnames2{i}) > 0
            fprintf(1,'Slave orbit file for pair %d is named %s\n',i,ofnames2{i});
        else
            error(sprintf('Cannot find  %80s\n',ofnames2{i}));
        end
    end
end

if nk == numel(pfnames)
    disp('Number of interferogram files (pairs) to adjust');
    np = length(imast)
else
    error(sprintf('Missing %d of %d interferogram files\n'...
        ,numel(pfnames)-nk,numel(pfnames)));
end

% approximate baseline in meters from Ha in meters
% valid for ERS only!!!
if abs(mpercy - 0.028) > 0.001
    warning('Approximating bperp\n');
end
bperp = 10000 ./ hamb;

% find the species
dispflag = 0;
[species, DD, tepochs, iepochs, iuniqorbs, uniqdates] = findspecies(t0,t1,dispflag,imast,mdate,islav,sdate);

disp('Number of distinct epochs overall'); me = length(tepochs)
disp('Number of species overall');         mfam = me - rank(DD)

[ns,ms] = size(species);
if ns ~= mfam
    error(sprintf('Error counting species: %d %d\n',ns,mfam));
end

% loop over species to count the number of pairs needed
np4 = 0;
for i=1:ns
    fprintf(1,        'Species %1s:',char(i+64));
    fprintf(fidtxtout,'Species %1s:',char(i+64));
    for j=1:ms
        if isfinite(species(i,j)) == 1
            fprintf(1,         ' %5d',species(i,j));
            fprintf(fidtxtout, ' %5d',species(i,j));
            if j > 1
                np4 = np4 +1;
            end
        end
    end
    fprintf(1,        '\n');
    fprintf(fidtxtout,'\n');
end

%    % loop over epochs to list pairs
np2 = 0;np3=0;
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
    
    [is,jm] = find(species == k);
    % is is index to species
    % jm is index of member within species
    for kk=1:numel(kp)
        np3 = np3+1;
        %imast(kp(kk))
        
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

% estimate the pseudoabsolute baselines
bpest = adjustbp(tepochs,DD,bperp, species,iuniqorbs, uniqdates);
nf=nf+1;h(nf)=figure;
plotbp(tepochs, bpest, DD, species, iuniqorbs, uniqdates, 0, 'orbital separation (Bperp) [m]');
%ktour = plotbp(tepochs, bpest, DD, species,iuniqorbs, uniqdates, 3);
%feval(printfun,sprintf('%s_%02d',mfilename,nf));
%feval(printfun,sprintf('%s_SPECIES',runname));
feval(printfun,sprintf('%s_SPECIES',runname));

if np > np4
    warning(sprintf('Number of pairs (%d) is larger than necessary (%d). Check for redundant pairs in list above.\n',np,np4));
end

%should estimate the pseuodabsolute doppler values
dops = zeros(size(tepochs));

% find out about the DEM
[isgeo,y1,x1,nl,nc,l1,c1,ml,mc,dl,dc,fi2,lat0,lon0,y0,x0,hemisphere,iutmzone] = read_dem_descriptor(demdescfile);

% Handle longitude west of Greenwich > 180
% if isgeo == 1
%    if xcenter > 180
%       xcenter = xcenter - 360;
%    end
%    if x1 > 180
%       x1 = x1 - 360;
%    end
% end

disp('number of columns in each interferogram');ncols = mc
disp('number of lines in each interferogram '); nrows = ml

if pselect == 5 || pselect == 7
    npix = ncols * nrows;
end
fprintf(1, 'Maximum number of pixels to be selected for inversion = %d\n',npix);

fprintf(1,'%12s %12s %12s\n',   'Location      ','X','Y');
fprintf(1,'%12s %12.0f %12.0f\n','DEM pixel 1  ',x1,y1);
% fprintf(1,'%12s %12.0f %12.0f\n','DEM extract 1',x1+dc*c1,y1+dl*l1);
% fprintf(1,'%12s %12.0f %12.0f\n','DEM extract N',x1+dc*(c1+mc),y1+dl*(l1+ml));
fprintf(1,'%12s %12.0f %12.0f\n','DEM extract 1',x1+dc*(c1-1)     ,y1+dl*(l1-1));
%fprintf(1,'%12s %12.0f %12.0f\n','DEM extract N',x1+dc*(c1-1+mc-1),y1+dl*(l1-1+ml-1));
% 2012-JAN-12
fprintf(1,'%12s %12.0f %12.0f\n','DEM extract N',x1+dc*(c1-1+mc),y1+dl*(l1-1+ml));
fprintf(1,'%12s %12.0f %12.0f\n','Center       ',xcenter, ycenter);

% find the pixel indices of the center of the subregion
% Corrected 2008-JUL-08    %% 26mar09lap [xy]center now in meters %%
% icenter = round((1e3*ycenter-y1)/dl - l1) % LINES   ARE COUNTED FROM FIRST LINE OF IMAGE
% jcenter = round((1e3*xcenter-x1)/dc - c1) % COLUMNS ARE COUNTED FROM LEFT (WEST) EDGE OF IMAGE
% 2012-JAN-11 next 2 lines are suspect
% icenter = round((ycenter-y1)/dl) - l1 + 1; % LINES   ARE COUNTED FROM FIRST LINE OF IMAGE
% jcenter = round((xcenter-x1)/dc) - c1 + 1; % COLUMNS ARE COUNTED FROM LEFT (WEST) EDGE OF IMAGE
icenter = round((ycenter-y1)/dl - l1); % LINES   ARE COUNTED FROM FIRST LINE OF IMAGE
jcenter = round((xcenter-x1)/dc - c1); % COLUMNS ARE COUNTED FROM LEFT (WEST) EDGE OF IMAGE
% icenter = ceil((ycenter-y1)/dl - l1); % LINES   ARE COUNTED FROM FIRST LINE OF IMAGE
% jcenter = ceil((xcenter-x1)/dc - c1); % COLUMNS ARE COUNTED FROM LEFT (WEST) EDGE OF IMAGE
% Handle case of simulated image
% if icenter - halfheight == 0
%     icenter = icenter+1;
% end
% if icenter + halfheight == ml
%     icenter = icenter-1;
% end
% if jcenter - halfwidth == 0
%     jcenter = jcenter+1;
% end
% if jcenter + halfwidth == mc
%     jcenter = jcenter-1;
% end
if icenter - halfheight < 1
    icenter = halfheight+1;
elseif icenter + halfheight > ml
    icenter = halfheight-1;
end
if jcenter - halfwidth < 1
    jcenter = halfwidth+1;
elseif jcenter + halfwidth > mc
    jcenter = halfwidth-1;
end

fprintf(1,'icenter = %d jcenter = %d\n',icenter,jcenter);

% % SANITY CHECK
% if icenter < 1 || icenter > mc
%     error(sprintf ('ERROR: icenter %d is out of bounds.\n',icenter));
% end
% if jcenter < 1 || jcenter > ml
%     error(sprintf ('ERROR: jcenter %d is out of bounds.\n',jcenter));
% end
% if icenter - halfheight < 1
%     error(sprintf ('ERROR: icenter (%d) minus halfheight (%d) is less than 1.\n',icenter,halfheight));
% end
% if icenter + halfheight > ml
%     error(sprintf ('ERROR: icenter (%d) plus halfheight (%d) is greater than ML (%d).\n',icenter,halfheight,ml));
% end
% if jcenter - halfwidth < 1
%     error(sprintf ('ERROR: jcenter (%d) minus halfwidth (%d) is less than 1.\n',jcenter,halfwidth));
% end
% if jcenter + halfwidth > mc
%     error(sprintf ('ERROR: jcenter (%d) plus halfwidth (%d) is greater than MC (%d).\n',jcenter,halfwidth,mc));
% end


% Define pixel indices of the sub-region here as a rectangle
% halfwidth and halfheight
% isub = (icenter-halfheight:icenter+halfheight) + 1;
% jsub = (jcenter-halfwidth :jcenter+halfwidth)  + 1;
% 2012-JAN-12 next 2 lines imply that subregion will have side lengths that
% are ODD numbers of pixels
if dl > 0
    isub = (icenter+halfheight:-1:icenter-halfheight);
    jsub = (jcenter-halfwidth :+1:jcenter+halfwidth);
else
    isub = (icenter-halfheight:+1:icenter+halfheight);
    jsub = (jcenter-halfwidth :+1:jcenter+halfwidth);
end
if dc < 0
    error('ERROR: NEGATIVE DC! PANIC STOP!');
end

% 2011-JUN-17 Kurt
% if min(isub) < 1 || max(isub) > nrows
%     isub = (icenter-halfheight:icenter+halfheight);
% end
% iok = intersect(find(isub > 0),find(isub <= nrows));
% if numel(iok) ~= 2*halfheight + 1
%     warning('Truncating isub.');
% end
% isub=isub(iok);
if min(isub) < 1
    warning(sprintf('ERROR: isub is less than 1 Check ycenter and halfheight.\n'));
    warning('Truncating isub LT 1');
    %   isub(isub < 1) = 1;
    isub(isub < 1) = NaN;
    nerrors = nerrors + 1;
end
if max(isub) > nrows
    warning(sprintf('ERROR: isub is greater than the number of rows in DEM. Check ycenter and halfheight.\n'));
    warning(sprintf('Truncating isub GT %d\n',nrows));
    %    isub(isub > nrows) = nrows;
    isub(isub > nrows) = NaN;
    nerrors = nerrors + 1;
end
isub = isub(find(isfinite(isub)));
isub = unique(isub);
%iok = intersect((xax >= x1+dc*(c1-1)),(xax <= x1+dc*(c1-1+mc)));
%isub = iok;
% if min(jsub) < 1 || max(jsub) > ncols
%     jsub = (jcenter-halfwidth :jcenter+halfwidth);
% end
% 2011-JUL-12
% jok = intersect(find(jsub > 0),find(jsub <= ncols));
% if numel(jok) ~= 2*halfwidth + 1
%     warning('Truncating jsub.');
% end
% jsub=jsub(jok);
%jok = intersect((yax >= y1+dl*(l1-1)),(yax <= y1+dl*(l1-1+ml)));
%jsub = jok;
if min(jsub) < 1
    warning(sprintf('ERROR: jsub is less than 1 Check xcenter and halfwidth.\n'));
    warning(sprintf('Truncating jsub LT %d\n',1));
    %jsub(jsub < 1) = 1;
    jsub(jsub < 1) = NaN;
    nerrors = nerrors + 1;
end
if max(jsub) > ncols
    warning(sprintf('ERROR: jsub is greater than number of columns in DEM. Check xcenter and halfwidth.\n'));
    warning(sprintf('Truncating jsub GT %d\n',ncols));
    %jsub(jsub > ncols) = ncols;
    jsub(jsub > ncols) = NaN;
    nerrors = nerrors + 1;
end
jsub = jsub(find(isfinite(jsub)));
jsub = unique(jsub);

nrsub = numel(isub);
ncsub = numel(jsub);
fprintf(1,'Subregion contains NRSUB = %d rows and NCSUB = %d columns\n',nrsub,ncsub);

% location for E-W profile
iprof = floor(halfheight); % index in sub-region
% location for NS profile
jprof = floor(halfwidth); % index in sub-region


% Easting and Northing in Meters
if isgeo >= 2  % CARTOGRAPHIC coordinates
    %     xax = (x1 + c1*dc + dc*[1:ncols])/1.0e3; % fencepost problem here?
    %     yax = (y1 + l1*dl + dl*[1:nrows])/1.0e3;
    %     xax = (x1 + dc*((c1-1) + (0:ncols) + 0.0)); % pixel registration
    %     yax = (y1 + dl*((l1-1) + (0:nrows) + 0.0));
    % 2012-JAN-12 Next 2 lines are correct
    %   xax = x1 + dc*((c1-1) + (0:ncols-1)); % pixel registration
    %  yax = y1 + dl*((l1-1) + (0:nrows-1));
    % 2012-JAN-14
    xax = x1 + dc*(c1 + (0:ncols-1)); % pixel registration
    yax = y1 + dl*(l1 + (0:nrows-1));
    if numel(xax) ~= ncols || numel(yax) ~= nrows
        warning(sprintf('Problem with axes: numel(xax) ncols  numel(yax) nrows %d %d %d %d\n'...
            ,numel(xax),ncols,numel(yax),nrows));
    end
    xstart = xax(1);
    ystart = yax(1);
    dx = dc;
    dy = dl;
    xcenter1 = xcenter;
    xcenter2 = (min(xax(jsub))+max(xax(jsub)))/2.0;
    xcenter3 = xax(jcenter);
    ycenter1 = ycenter;
    ycenter2 = (min(yax(isub))+max(yax(isub)))/2.0;
    ycenter3 = yax(icenter);
    
    
    % Get origin of UTM projection
    if isfinite(lat0) == 0 || isfinite(lon0) == 0
        %utmzone0 = '02 U'; % Okmok
        if abs(iutmzone) > 0 && numel(hemisphere) > 0
            if strcmpi(hemisphere,'N') == 1
                utmzone0 = sprintf('%02d X',iutmzone);
            elseif strcmpi(hemisphere,'S') == 1
                utmzone0 = sprintf('%02d X',iutmzone);
            else
                warning(sprintf('Unknown zone for UTM projection. Setting to 31 N\n'));
                utmzone0 = '31 N'; % just east of Greenwich, just north of equator
                nerrors = nerrors + 1;
            end
        else
            warning(sprintf('Unknown zone for UTM projection. Setting to 31 N\n'));
            utmzone0 = '31 N'; % just east of Greenwich, just north of equator
            nerrors = nerrors + 1;
        end
        [lat0, lon0] = utm2deg(xcenter,ycenter,utmzone0);
    end
    
    % Get origin of UTM projection in degrees
    % signed longitude
    if lon0 > 180
        % slon0 = lon0-180; % corrected 2010-04-25 Kurt
        slon0 = lon0-360;
    else
        slon0 = lon0;
    end
    [xutm0,yutm0,utmzone0] = deg2utm(lat0,slon0);
    
    switch isgeo
        case 2
            % Get origin of our private projection
            if isfinite(x0) == 0
                x0 = xutm0;
            end
            if isfinite(y0) == 0
                y0 = yutm0;
            end
            [latcenter, loncenter] = utm2deg(xcenter,ycenter,utmzone0);
            
            % approximate pixel dimensions in degrees
            [lattemp, lontemp  ] = utm2deg(xcenter+dx,ycenter,utmzone0);
            dlon = lontemp-loncenter  % step size in longitude in degrees
            [lattemp, lontemp  ] = utm2deg(xcenter    ,ycenter+dy,utmzone0);
            dlat = lattemp - latcenter  % step size in latitude in degrees (can be negative)
        case 3
            fprintf(1,'Latitude of  Lambert origin %20.10f\n',lat0);
            fprintf(1,'Longitude of Lambert origin %20.10f\n',slon0);
            fprintf(1,'UTM zone of Lambert origin %s\n',utmzone0);
            latcenter = lat0;
            loncenter = lon0;
            warning(sprintf('ROUGHLY approximating pixel dimensions in degrees\n'));
            dlat = 2.0*pi*dl*6371.0E3/360.0  % step size in latitude in degrees (can be negative)
            dlon = 2.0*pi*dc*6371.0E3/360.0*cos(pi*latcenter/180.)  % step size in longitude in degrees (can be negative)
        otherwise
            error(sprintf('Unknown isgeo %d\n',isgeo));
    end
else           % GEOGRAPHIC (longitude, latitude) coordinates in degrees
    %    lonax = (x1 + c1*dc + dc*[1:ncols]); % fencepost problem here?
    %    latax = (y1 + l1*dl + dl*[1:nrows]);
    %    lonax = sort(lonax);
    %    latax = sort(latax);
    %     lonax = (x1 + c1*dc + dc*[1:ncols+1]); %
    %     latax = (y1 + l1*dl + dl*[1:nrows+1]);
    if dl > 0 % 20120924 for Chile with Helene's DEM
        lonax = x1 + dc*((c1-1) + (0      :+1:ncols-1) + 0.0); % pixel registration
        latax = y1 + dl*((l1-1) + (nrows-1:-1:0      ) + 0.0); % pixel registration
    else % usual case with decreasing latitude
        lonax = x1 + dc*((c1-1) + (0:ncols-1) + 0.0); % pixel registration
        latax = y1 + dl*((l1-1) + (0:nrows-1) + 0.0); % pixel registration
    end
    if dc < 0
        error('ERROR: NEGATIVE DC! PANIC STOP!');
    end
    %    lonax = linspace(x1+(c1-1)*dc,x1+mc*dc,ncols);
    %    latax = linspace(y1+(l1-1)*dl,y1+ml*dl,nrows);
    
    % Old UTM projection
    %    [xstart, ystart,utmzone1] = wgs2utm(latax(1)  ,lonax(1));
    %    [xend,   yend,  utmzone2] = wgs2utm(latax(end),lonax(end));
    % 2009-JUL-20 use new UTM projection
    [xstart, ystart,utmzone1] = deg2utm(latax(1)  ,lonax(1));
    %  [xend,   yend,  utmzone2] = deg2utm(latax(end),lonax(end));
    [xend,   yend,  utmzone2] = deg2utm(latax(1)+dl,lonax(1)+dc);
    dx = xend - xstart % East-West   dimension of pixel in meters
    dy = yend - ystart % North-south dimension of pixel in meters
    
    dlon = dc;  % step size in longitude in degrees
    dlat = dl;  % step size in latitude in degrees (can be negative)
    %     xax = (xstart + dx*[1:ncols+1]);
    %     yax = (ystart + dy*[1:nrows+1]);
    xax = xstart + dx*(0:ncols-1); % pixel registration
    yax = ystart + dy*(0:nrows-1);
    
    %    [xstart, ystart,utmzone1] = deg2utm(lonax(1)  ,latax(1));
    %    [xend,   yend,  utmzone2] = deg2utm(lonax(end),latax(end));
    if strncmp(utmzone1,utmzone2,4) == 0
        warning(sprintf('WARNING: UTM zones do not match for corners of subregion. %s %s\n'...
            ,utmzone1,utmzone2));
    end
    %    xax = linspace(xstart, xend, ncols+1);
    %    yax = linspace(ystart, yend, nrows+1);
    %    xax = sort(xax);
    %    yax = sort(yax);
    %    xax = xax(1:end-1) + (xax(2)-xax(1))/2;
    %    yax = yax(1:end-1) + (yax(2)-yax(1))/2;
    %     xax = xax/1.0e3;
    %     yax = yax/1.0e3;
    %
    % Convert (lat, lon) to meters (not km). TRICKY! Kurt 2009-APR-01
    %[xcenter1, ycenter1,utmzone0] = wgs2utm(ycenter ,xcenter) Kurt
    %2010-JUL-07    
    loncenter = x1 + dc * (jcenter-1);
    latcenter = y1 + dl * (icenter-1);
    [xcenter1, ycenter1,utmzone0] = deg2utm(latcenter ,loncenter)
    
    if numel(xcenter1) < 1 || numel(ycenter1) < 1
        latcenter
        loncenter
        xcenter1
        ycenter1    
        error('xcenter1 or ycenter1 is not defined');
    end
    
    % Check by approximation
    xcenter2 = xstart + dx*(xcenter-(x1 + (c1-1)*dc))/dc
    ycenter2 = ystart + dy*(ycenter-(y1 + (l1-1)*dl))/dl
    xcenter3 = xstart + dx * (jcenter - 1)
    ycenter3 = ystart + dy * (icenter - 1)
    
    if abs(xcenter2 - xcenter1) < dx && abs(ycenter2 - ycenter1) < dy
        xcenter = xcenter1;
        ycenter = ycenter1;
    else
        warning('WARNING: problem with UTM coordinates for center point.');
        format bank
        xcenter1
        xcenter2
        diffx=xcenter2-xcenter1
        ycenter2
        ycenter1
        diffy=ycenter2-ycenter1
        format compact
        xcenter = xcenter2
        ycenter = ycenter2
        %nerrors = nerrors + 1;
    end
    
    if strcmp(utmzone0,utmzone1) == 0
        warning(sprintf('WARNING: UTM zones do not match for center of subregion. %s %s\nApproximating...\n'...
            ,utmzone0,utmzone1));
        format bank
        xcenter = xcenter2
        ycenter = ycenter2
        format compact
        %nerrors = nerrors + 1;
    end
    
end


% extrema of extracted area
xmin = min(xax);
xmax = max(xax);
ymin = min(yax);
ymax = max(yax);

% READ THE DEM
i2dem = read_i2(fi2,nc);
[nrdem,ncdem]=size(i2dem);
fprintf(1,'I2DEM has %d rows and %d columns\n',nrdem,ncdem);
figure;imagesc(i2dem);
if nrdem >= ml && ncdem >= mc
    % fill with zeros
    zframe=zeros(ml,mc);
    % use values from DEM where available
    izl1 = l1;
    izc1 = c1;
    izl2 = min([(l1+ml-1),nrdem]);
    izc2 = min([(c1+mc-1),ncdem]);
    jzl1 = 1;
    jzl2 = izl2-izl1+1
    jzc1 = 1;
    jzc2 = izc2-izc1+1
    zframe(jzl1:jzl2,jzc1:jzc2)=double(i2dem(izl1:izl2,izc1:izc2));
    %try to flip dem
    if dl>0
        warning('DL > 0 : flipping zframe\n');
        zframe=flipud(zframe);
    end
else
    zframe = double(i2dem(1:nrdem,1:ncdem));
    warning(sprintf('WARNING: DEM wrong size. Expected %d rows and %d columns. Setting undefined elevations to 0 meters.\n',ml,mc));
    %zframe(isub,jsub)=zeros(numel(isub),numel(jsub));
end
disp('Nrows, Ncols of Extracted Region in DEM.');size(zframe)

% Decide how to select pixels
switch pselect
    case 0
        fprintf(1,'Selecting all pixels.\n');
    case 1
        fprintf(1,'Selecting pixels randomly.\n');
    case 2
        fprintf(1,'Reading previous pixel indices from existing files ikeep.mat and jkeep.mat.\n');
    case 3
        fprintf(1,'Selecting pixels using quadtree subsampling, then randomly\n');
    case 4
        fprintf(1,'Selecting pixels based on coherence.\n');
    case 5
        fprintf(1,'Selecting pixels using pha2qls resampling.\n');
    case 6
        fprintf(1,'Selecting pixels based on quadtree of stacked phase.\n');
    case 7
        fprintf(1,'Selecting pixels using pha2qls. Values are range gradient.\n');
    case 9
        fprintf(1,'Selecting pixels randomly. Values are range gradient.\n');
    otherwise,
        error(sprintf('Unkown value of pselect (%d)\n',pselect));
end
switch pselect
    case 0
        npix = nrsub*ncsub;
        fprintf(1,'Selecting all %d pixels.\n',npix);
        [ikeep,jkeep] = meshgrid(icenter-halfheight:icenter+halfheight,jcenter-halfwidth:jcenter+halfwidth);
        ikeep = colvec(ikeep);
        jkeep = colvec(jkeep);
        disp 'ikeep';size(ikeep)
        disp 'jkeep';size(jkeep)
    case {1,3} % Choose pixels randomly
        if npix > nrsub * ncsub
            warning(sprintf('NPIX (%d) GT NRSUB (%d) * NCSUB (%d) = %d\nResetting\n'...
                ,npix,nrsub,ncsub,nrsub*nsub));
            npix =  nrsub * ncsub
        end
        fprintf(1,'Selecting %d pixels randomly\n',npix);
        %20130624 - fprintf(1,'Initializing random number generator.\n');
        %        Replace the default stream with a stream whose seed is based on CLOCK, so
        %        RAND will return different values in different MATLAB sessions.  NOTE: It
        %        is usually not desirable to do this more than once per MATLAB session.
        % 20130624 - Do something about the error message below.
        %         RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
        %         Warning: The RandStream.setDefaultStream static method will be removed in a future
        %         release.  Use RandStream.setGlobalStream instead.
        %         > In RandStream.RandStream>RandStream.setDefaultStream at 456
        %         In gipht_step1 at 691
        %         In gipht at 119
        %20130624 - My understanding is that if we do NOT initialize, then the RNG
        %will generate a different number in each session of GIPHT.


        %         ikeep = icenter + floor(2*halfheight*(rand([1,npix])-0.5)) + 1;
        %         jkeep = jcenter + floor(2*halfwidth *(rand([1,npix])-0.5)) + 1;
        % 2012-JAN-12
        %             ikeep = icenter + ceil((2*halfheight+1)*(rand([1,npix])-0.5));
        %             jkeep = jcenter + ceil((2*halfwidth +1)*(rand([1,npix])-0.5));
        irand = halfheight*2*(rand([1,npix])-0.5);
        jrand = halfwidth *2*(rand([1,npix])-0.5);
        ikeep = icenter + floor(abs(irand)) .* sign(irand);
        jkeep = jcenter + floor(abs(jrand)) .* sign(jrand);
        fprintf(1,'Saving %d pixel indices to ikeep.mat and jkeep.mat.\n',npix);
        save ikeep ikeep; save jkeep jkeep
    case 2
        % See if files exist
        fdi = fopen ('ikeep.mat','r'); fdj = fopen ('jkeep.mat','r');
        if fdi == -1 || fdj == -1
            error(sprintf('Could not open ikeep.mat and/or jkeep.mat\nPlease set select = 1.\n'));
        else
            fclose(fdi);fclose(fdj);
            fprintf(1,'Reading previous pixel indices from existing files ikeep.mat and jkeep.mat.\n');
            load ikeep; load jkeep
            if length(ikeep) ~= npix || length(jkeep) ~= npix
                error(sprintf(1, ...
                    'Number of pixels (%d) does not match number in [ij]keep.mat.\nPlease set select = 1.\n',npix));
            end
        end
    case {4,5,6,7,9}
        ikeep(1) = icenter - halfheight;
        jkeep(1) = jcenter - halfwidth;
        ikeep(2) = icenter + halfheight;
        jkeep(2) = jcenter + halfwidth;
    otherwise,
        error(sprintf('Unkown value of pselect (%d)\n',pselect));
end



fprintf(1,'Extrema of coordinates of images (extracted area)\n');
fprintf(1,'%12s %12s %12s\n','Location','X','Y');
fprintf(1,'%12s %12.0f %12.0f\n','Min '   ,xmin,ymin);
fprintf(1,'%12s %12.0f %12.0f\n','Max'    ,xmax,ymax);
fprintf(1,'%12s %12.0f %12.0f\n','Diff '  ,xmax-xmin, ymax-ymin);
fprintf(1,'%12s %12.0f %12.0f\n','Center' ,xcenter,ycenter);

xsubmin = min(xax(jsub));
xsubmax = max(xax(jsub));
ysubmin = min(yax(isub));
ysubmax = max(yax(isub));

fprintf(1,'Extrema of coordinates of subregion\n');
fprintf(1,'%12s %12s %12s\n','Location','Easting(m)','Northing(m)');
fprintf(1,'%12s %12.0f %12.0f\n','Min    ',xsubmin,         ysubmin);
fprintf(1,'%12s %12.0f %12.0f\n','Max    ',xsubmax,         ysubmax);
fprintf(1,'%12s %12.0f %12.0f\n','Diff   ',xsubmax-xsubmin, ysubmax-ysubmin);
fprintf(1,'%12s %12.0f %12.0f\n','Center ',xcenter,         ycenter);
fprintf(1,'%12s %12.0f %12.0f\n','Center1',xcenter1,        ycenter1);
fprintf(1,'%12s %12.0f %12.0f\n','Center2',xcenter2,        ycenter2);
fprintf(1,'%12s %12.0f %12.0f\n','Center3',xcenter3,        ycenter3);
fprintf(1,'%12s %12.0f %12.0f\n','Max-Cen ',xsubmax-xcenter, ysubmax-ycenter);
fprintf(1,'%12s %12.0f %12.0f\n','Cen -Min',xcenter-xsubmin, ycenter-ysubmin);
fprintf(1,'%12s %12.0f %12.0f\n','Max-Cen3',xsubmax-xcenter3, ysubmax-ycenter3);
fprintf(1,'%12s %12.0f %12.0f\n','Cen3-Min',xcenter3-xsubmin, ycenter3-ysubmin);

if xcenter < xsubmin || xcenter > xsubmax
    warning('ERROR xcenter outside subregion.');
    nerrors = nerrors + 1;
end
if ycenter < ysubmin || ycenter > ysubmax
    warning('ERROR ycenter outside subregion.');
    nerrors = nerrors + 1;
end

% plot subregion on DEM
nf=nf+1;h(nf)=figure;
imagesc([xmin, xmax]/1e3,[ymax ymin]/1e3,double(zframe));
hold on;axis xy; axis tight;
plot([xsubmin xsubmax xsubmax xsubmin xsubmin]/1e3...
    ,[ysubmin ysubmin ysubmax ysubmax ysubmin]/1e3...
    ,'k+-');
xlabel('Easting (km)');ylabel('Northing (km)');
title(sprintf('DEM - elevation in meters %s %s NL = %d NC = %d ML = %d MC = %d'...
    ,demdescfile,fi2 ...
    ,nl,nc,ml,mc)...
    ,'Interpreter','None');
colorbar;cmapblacknan;
feval(printfun,sprintf('dem'));

% make new descriptor file
% starting line, column indices for subregion
% 2012-JAN-11 The next two lines are suspect
% l1_sub = l1 - 1 + icenter - halfheight;
% c1_sub = c1 - 1 + jcenter - halfwidth;
% 2012-JAN-11
% l1_sub = l1  + icenter - halfheight;
% c1_sub = c1  + jcenter - halfwidth;
% 2012-JAN-11
% l1_sub = l1 - 1 + icenter -1 - halfheight;
% c1_sub = c1 - 1 + jcenter -1 - halfwidth;
% % This will only work for isgeo > 2
if isgeo > 2
    if dl > 0
        l1_sub = round(l1-1 + (ysubmin - y0)/dl);
    else
        l1_sub = round(l1-1 + (ysubmax - y0)/abs(dl));
    end
    if dc > 0
        c1_sub = round(c1-1 + (xsubmin - x0)/dc);
    else
        c1_sub = round(c1-1 + (xsubmax - x0)/abs(dc));
    end
else
    l1_sub = l1 - 1 + icenter - halfheight;
    c1_sub = c1 - 1 + jcenter - halfwidth;   
end
% if dl > 0
%     l1_sub = l1 - 1 + find(abs(yax-ysubmin) < abs(dl));
% else
%     l1_sub = l1 - 1 + find(abs(yax-ysubmax) < abs(dl));
% end
% c1_sub = c1 - 1 + find(abs(xax-xsubmin) < abs(dc));


% number of lines, columns in subregion 
% 20130624make symmetry work according to odd or even number
if mod(nrsub,2) == 1
    ml_sub = 2*halfheight+1;
else
    ml_sub = 2*halfheight+2;
end
if mod(ncsub,2) == 1
    mc_sub = 2*halfwidth+1;
else
    mc_sub = 2*halfwidth+2;
end
% 2012-06-24 even number
% ml_sub = 2*halfheight;
% mc_sub = 2*halfwidth;

% if abs((1+max(ikeep)-min(ikeep)) - ml_sub) > 0
%     error(sprintf('Height of ikeep = %d DIFFERS FROM Height of subregion = %d\n',1+max(ikeep)-min(ikeep),ml_sub));
% end
% if abs((1+max(jkeep)-min(jkeep)) - mc_sub) > 0
%     error(sprintf('Width  of jkeep = %d DIFFERS FROM Width  of subregion = %d\n',1+max(jkeep)-min(jkeep),mc_sub));
% end
if abs(nrsub - ml_sub) > 0
    warning(sprintf('nrsub = %d DIFFERS FROM Height of subregion = %d\n',nrsub,ml_sub));
    nerrors = nerrors + 1;
end
if abs(ncsub - mc_sub) > 0
    warning(sprintf('ncsub = %d DIFFERS FROM Width  of subregion = %d\n',ncsub,mc_sub));
    nerrors = nerrors + 1;
end

% copy the input descriptor file to the output descriptor file, changing only the parameters above
ierr = write_dem_descriptor(demdescfile,demdescfile2,isgeo,y1,x1,nl,nc,l1_sub,c1_sub,ml_sub,mc_sub,dl,dc,fi2,lat0,lon0,y0,x0,hemisphere,iutmzone);

if ierr == 0
    fprintf(1,'Wrote DEM descriptor for subregion to file named %s\n',demdescfile2);
else
    error(sprintf('Could not write DEM descriptor for subregion to file named %s\n',demdescfile2));
end

% make a stack
if  pselect == 6
    sphnam =sprintf('%s_pstackpha.pha',runname); % output phase stack
    qphnam =sprintf('%s_qstackpha.pha',runname); % quadtreed phase stack
    qlsnam =sprintf('%s_qstack_pha.i2',runname); % quadtree indices
    
    % perform stacking
    np0 = stack_pha(pfnames,ncols,nrows,sphnam);
    
    ierr2 = write_pst(fnamein,fitfun,mparam,p0,p1,psig,pnames,bounds);
    sph=read_pha(sphnam,ncols);
    figure;imagesc(sph);colorbar;cmapblackzero;
    title('wrapped phase stack (256 DN per cycle)')
    xlabel('column index');ylabel('row index');
    
    % perform quad-tree partitioning on the stack
    ithreshS = 32; % threshold circular mean deviation in DN [0 127]
    pixinpatchS = 9;   % mininum number of good pixels in a quad
    nquad1 = pha2qls(sphnam,ncols,nrows,qphnam,qlsnam,ithreshS,pixinpatchS);
    
    qph=read_pha(qphnam,ncols);
    figure;imagesc(qph);colorbar;cmapblackzero;
    title('wrapped phase stack after quad-tree partitioning (256 DN per cycle)')
    xlabel('column index');ylabel('row index');
end

if (pselect == 2 || fnewer(orbsavefile,'gipht.in') == 0) && fexist(orbsavefile) == 1
    if pselect == 2 && fexist(orbsavefile) == 1
        fprintf(1,'Using existing file named %s containing orbit information...\n',orbsavefile);
        load(orbsavefile);
        orbits_loaded = 1;
    end
end

% set pointers for indices of pixels
ippix1 = zeros(np,1); i1 = 0; % index to first pixel in pair
ippix2 = zeros(np,1); i2 = 0; % index to last pixel in pair

% loop over pairs
nk=0;ndata1=0;i1=0;i2=0;kk=0;

for i = 1:np
    npixinpair=0;
    
    % index master and slave
    kmast = find(DD(i,:) == -1);
    kslav = find(DD(i,:) == +1);
    
    %read phase data from interferogram
    %2009-JUN-18 phaimg = read_pha(fn0,ncols)/256; % in cycles
    
    fn0 = pfnames{i};
    nbytes = fsize(fn0);
    if  nbytes > 0
        phaimg = read_pha(fn0,ncols,nrows); % in DN [-128, 127]
    else
        phaimg = [];
        warning(sprintf('Phase file named %s is non-existant or empty\n',fn0));
        nerrors = nerrors + 1;
    end
    
    %figure;imagesc(phaimg);axis ij;xlabel('column index J');ylabel('row index I');title('before flipud');hold on;plot(jcenter,icenter,'k*');
    %helene try to flip the image
    if dl>0
        phaimg=flipud(phaimg);
        warning('DL > 0 flipping phaimg');
        %figure;imagesc(phaimg);axis ij;xlabel('column index J');ylabel('row index I');title('after flipud');hold on;plot(jcenter,icenter,'k*');        
    end
    nbytes = numel(phaimg);
    if nbytes == nrows*ncols
        nk=nk+1; % count it as a success
    else
        warning(sprintf('Number of bytes (%d) in phase file %s is not equal to product of nrows (%d) and ncols (%d) = (%d)'...
            ,nbytes,fn0,nrows,ncols,nrows*ncols));
        nerrors = nerrors + 1;
    end
    
    % read the coherence file if needed
    if pselect == 4
        fnc = strrep(strrep(fn0,'smp','coh'),'.pha','.byt');
        cohimg = read_oct(fnc,ncols);
    end
    
    % simulate random noise
    % phar = 0.1+0.2*randn([nrows,ncols]);
    % phaimg = angle(complex(cos(phar*2*pi),sin(phar*2*pi)))/2/pi;
    
    titlestr = sprintf('Phase: Pair %3d %s orbs %5d %5d years %7.1f to %7.1f Dt = %.4f yr NROWS = %d by NCOLS= %d'...
        ,i,strrep(fn0,'_','\_'),iuniqorbs(kmast),iuniqorbs(kslav)...
        ,tepochs(kmast),tepochs(kslav),tepochs(kslav)-tepochs(kmast)...
        ,nrows,ncols);
    
    % display the whole image
    nf=nf+1;h(nf)=figure;
    
%     if dl < 0
%         imagesc([xmin, xmax]/1e3,[ymax ymin]/1e3,double(phaimg)/256,[-0.5 0.5]);
%     else
%         warning('DL > 0: reversing order of ymin and ymax');
%         imagesc([xmin, xmax]/1e3,[ymin ymax]/1e3,double(phaimg)/256,[-0.5 0.5]);
%     end
    imagesc([xmin, xmax]/1e3,[ymax ymin]/1e3,double(phaimg)/256,[-0.5 0.5]);
    ctab = cmapblackzero(1); % black at zero at bottom of color bar
    colormap(ctab);
    hold on;
    plot([xsubmin xsubmax xsubmax xsubmin xsubmin]/1e3...
        ,[ysubmin ysubmin ysubmax ysubmax ysubmin]/1e3,'k+-');
    xlabel('Easting (km)');ylabel('Northing (km)');
    axis xy; axis tight; axis equal;hold on;cmapblackzero;
    colorbar; title(titlestr);
    %feval(printfun,sprintf('%s_%02d',mfilename,nf));
    
    if nerrors > 0
        error(sprintf('ERROR: too many errors (%d) to continue!\n',nerrors));
    end
    
    if ismember(pselect,[3,5,6,7,9])
        % perform quadtree partitioning using pha2qls program
        % WITH GRADIENTS
        qlsnam=strrep(fn0,'.pha','.qls');
        if numel(strfind(fn0,'psp_')) > 0
            qspnam=strrep(fn0,'psp_','qsp_');
            grxnam=strrep(fn0,'psp_','grx_');grxnam=strrep(grxnam,'.pha','.i2');
            grynam=strrep(fn0,'psp_','gry_');grynam=strrep(grynam,'.pha','.i2');
        else
            qspnam='qsp.pha';
            grxnam='grx.i2';
            grynam='gry.i2';
        end
        
        qspnames{i}=qspnam;
        grxnames{i}=grxnam;
        grynames{i}=grynam;
        nqcols = 6;
        
        % perform quadtree resampling
        if fexist(qspnam) <= 0 || fexist(qlsnam) <= 0 || fexist(grxnam) <= 0 || fexist(grynam) <= 0 ...
                || fnewer(qspnam,'gipht.in') == 1 ...
                || fnewer(qspnam,fn0)     == 1 ...
                || fsize(qspnam) ~= fsize(fn0)
            nquad1 = pha2qls(fn0,ncols,nrows,qspnam,grxnam,grynam,qlsnam,ithresh,pixinpatch,maxcmd,maxpix,pha2qlsname);
        end
        
        % read binary files in DN
        %psp= read_pha(fn0,ncols);
        psp = phaimg; % already read above
        if size(psp) ~= [nrows,ncols]
            size(psp)
            error('PSP wrong size');
        end
        qsp= read_pha(qspnam,ncols,nrows);
        if size(qsp) ~= [nrows,ncols]
            size(qsp)
            error('QSP wrong size');
        end
        grx= read_i2(grxnam,ncols);
        gry= read_i2(grynam,ncols);
        if size(grx) ~= size(gry)
            size(grx)
            size(gry)
            error('GRX and GRY different sizes');
        end
        
        % transfer missing values (coded as zeros) from phase image
        inan=find(qsp==0); % test while still integers
        psp(inan)=NaN;
        qsp(inan)=NaN;
        grx(inan)=NaN;
        gry(inan)=NaN;
        
        % change from integers to doubles
        psp=2*pi*double(psp)/256.0;% convert to radians
        qsp=2*pi*double(qsp)/256.0;% convert to radians
        grx=2*pi*double(grx)/256.0/256.0; % convert to radians per pixel
        gry=2*pi*double(gry)/256.0/256.0; % convert to radians per pixel
        
        % show us what we have got
        figure;hold on;imagesc(psp/DNPC,[-0.5,0.5]);colorbar;
        ctab = cmapblackzero; % black at zero at middle of color bar
        colormap(ctab);
        if dl < 0; axis ij; else axis xy; end
        title('wrapped phase (cycles)');
        xlabel('column index');ylabel('row index');
        
        maxabs=0.5;
        figure;hold on;imagesc(qsp/DNPC,[-maxabs,maxabs]);colorbar;
        ctab = cmapblackzero; % black at zero at middle of color bar
        colormap(ctab);
        title('wrapped phase after quad-tree partitioning (cycles)')
        xlabel('column index');ylabel('row index');
        if dl < 0; axis ij; else axis xy; end
        
        if ismember(pselect,[7,9])
            maxabs = max(max(abs(grx/DNPC)));
            figure;imagesc(grx/DNPC,[-maxabs,maxabs]);colorbar;cmapblackzero;
            title('IMAGESC: East component of phase gradient after quad-tree partitioning (cycles/pixel)')
            xlabel('column index');ylabel('row index');
        end
    end
    
    if ismember(pselect,[5,6,7,9])
        % read the list of indices pointing to pixels in quads
        qlist = read_i2(qlsnam,nqcols);
        % first row is special
        % 20121004 if qlist(1,3) ~= ncols  || qlist(1,4) ~= nrows
        % first row is special because it contains number of rows and columns
        if typecast(qlist(1,3:4),'int32') ~= ncols  || typecast(qlist(1,5:6),'int32') ~= nrows
            error(sprintf('Number of columns (%d %d) or rows (%d %d) incorrect.\n'...
                ,typecast(qlist(1,3:4),'int32'),ncols...
                ,typecast(qlist(1,5:6),'int32'),nrows));
        end
        
        %        qi1=qlist(2:end,1);% column index of first pixel in quad
        %        qj1=qlist(2:end,2);% row index of first pixel in quad
%         %      THE TWO LINES ABOVE ARE WRONG! THE TWO LINES BELOW ARE CORRECT. 2010-JUL-08
%         qi1=qlist(2:end,2);% ROW index of first pixel in quad
%         qj1=qlist(2:end,1);% COL index of first pixel in quad
%         qkw=qlist(2:end,3);% width of square quad
%       20140106 Matlab needs these pointers to be real values
        qi1=double(qlist(2:end,1));  % Index to col of first pixel in patch
        qj1=double(qlist(2:end,2));  % Index to row of first pixel in patch
        qkw=double(qlist(2:end,3));% width of square quad

        %        if pselect == 5
        %            qqp=2*pi*double(qlist(2:end,4))/256;% phase value
        %        end
        %        if pselect == 7
        % 2011-MAR-24 - GREAT BIG BUG - NOW fixed in pha2qls3.c
        qqp=2*pi*double(qlist(2:end,4))/256/256;% phase value
        %      qqp=2*pi*double(qlist(2:end,4))/256;% phase value in RADIANS
        qgx=2*pi*double(qlist(2:end,5))/256/256;  % X-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
        qgy=2*pi*double(qlist(2:end,6))/256/256;  % Y-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
        %      end
        
        % only take patches larger than a single pixel
        %iok=find(qkw > 1);
        % only take patches with at least mininum
        iok=find((qkw.*qkw) >= pixinpatch);
        qi1=qi1(iok);
        qj1=qj1(iok);
        qqp=qqp(iok);
        qkw=qkw(iok);
        qgx=qgx(iok);
        qgy=qgy(iok);
        
        %        %% sort in decreasing order of width
        %        [qkw,isort]=sort(qkw,1,'descend');
        %        qi1=qi1(isort);
        %        qj1=qj1(isort);
        %        qqp=qqp(isort);
        %        qgx=qgx(isort);
        %        qgy=qgy(isort);
        
        %shuffle order randomly  - but create an error with unit vector
        if pselect == 9
            jj=colvec(1:numel(qqp));
            iok=ishuffle(jj);
            qi1=qi1(iok);
            qj1=qj1(iok);
            qqp=qqp(iok);
            qkw=qkw(iok);
            qgx=qgx(iok);
            qgy=qgy(iok);
        end
        
        % build other indices
        % 2012-JAN-14 next 2 lines are correct according to
        % demoQ/test_quadphase9.m
        qi2=qi1+qkw-1;   % column index of last pixel in patch
        qj2=qj1+qkw-1;   % row index of last pixel in patch
        %        qim=double((qi1+qi2)/2.0);  % column index of middle of patch
        %        qjm=double((qj1+qj2)/2.0);  % row    index of middle of patch
        %      2011-JUN-22
        qim=round(double(qi1+qi2)/2.0);  % column index of middle of patch
        qjm=round(double(qj1+qj2)/2.0);  % row    index of middle of patch
        % 2012-JAN-14
        %         qim=floor(double(qi1+qi2)/2.0);  % column index of middle of patch
        %         qjm=floor(double(qj1+qj2)/2.0);  % row    index of middle of patch
        %      2014-JAN-06
        %         qim=double(qi1+qi2)/2.0;  % column index of middle of patch
        %         qjm=double(qj1+qj2)/2.0;  % row    index of middle of patch

        
        % convert indices from C to Fortran convention
        qim = qim+1; qjm = qjm+1;
        qi1 = qi1+1; qi2 = qi2+1;
        qj1 = qj1+1; qj2 = qj2+1;
        
        % choose the observable
        switch pselect
            case 5
                qqq = qqp/DNPC;
                tstr='Phase after quad-tree partitioning';
                zstr='cycles'
            case  {7,9}
                qqq = qgx/DNPC;
                tstr='East component of phase gradient after quad-tree partitioning';
                zstr='Cycles/pixel'
            otherwise
                qqq = NaN;
                error(sprintf('Unknown value of pselect %d\n',pselect));
        end
        
        nquad=numel(qim);
        qskip = 1;
        ntake = npix;
        
        %        if nquad > npix
        %            fprintf(1,'WARNING: number of quadtree values (%d) exceeds number of pixels specified (%d). Taking first %d pixels.\n'...
        %                ,nquad,npix,ntake);
        %        end
        
        clear qim2 qjm2 qqp2
        
        k =0;
        while npixinpair < ntake && k < nquad - qskip
            k=k+qskip;
            %iphase = int8(0);
            rphase = 0; % double radians 2010-JAN-11
            
%             % 2012-JAN-19 MAJOR REMODELING WORK HERE
%             % Get Easting, Northing and Elevation coordinates in meters
%             if isgeo >= 2 % cartographic
%                 xxxm = x1 + dc * (c1-1 + double(qj1(k)-1) + double(1+qkw(k))/2.0);     % easting in meters
%                 yyym = y1 + dl * (l1-1 + double(qi1(k)-1) + double(1+qkw(k))/2.0);     % northing in meters
%             else
%                 xxxm =  xstart + dx * (double(qj1(k)-1) + double(1+qkw(k))/2.0);     % easting in meters
%                 yyym =  ystart + dy * (double(qi1(k)-1) + double(1+qkw(k))/2.0);     % northing in meters
%             end
            % 2014-JAN-06 Permute indices
            % Get Easting, Northing and Elevation coordinates in meters
            if isgeo >= 2 % cartographic
                xxxm = x1 + dc * (c1 + qi1(k) + qkw(k)/2.0);     % easting in meters
                yyym = y1 + dl * (l1 + qj1(k) + qkw(k)/2.0);     % northing in meters
            else
                xxxm =  xstart + dx * (qi1(k) + qkw(k)/2.0);     % easting in meters
                yyym =  ystart + dy * (qj1(k) + qkw(k)/2.0);     % northing in meters
            end
            
            %             % elevation in meters from DEM
            %             % 20130-05-13
            %             if qim(k) >= 1 && qim(k) <= nrdem && qjm(k) >= 2 && qjm(k) <= ncdem
            %                 zzzm = zframe(qim(k),qjm(k));
            %                 % topographic relief in east direction
            %                 % elevation difference
            %                 % in meters between adjacent pixels
            %                 if qjm(k) < mc
            %                     dzzzm = double(zframe(qim(k),qjm(k)+1)-zframe(qim(k),qjm(k)));
            %                 else
            %                     dzzzm = double(zframe(qim(k),qjm(k))-zframe(qim(k),qjm(k)-1));
            %                 end
            %             else
            %                 zzzm = 0;
            %                 dzzm = 0;
            %             end
            % elevation in meters from DEM
            % 2014-JAN-06 Permute indices
            % 20130-05-13
            if qim(k) >= 1 && qjm(k) <= nrdem && qim(k) >= 1 && qim(k) <= ncdem
                zzzm = zframe(qjm(k),qim(k));
                % topographic relief in east direction
                % elevation difference
                % in meters between adjacent pixels
                if qim(k) < mc
                    dzzzm = double(zframe(qjm(k),qim(k)+1)-zframe(qjm(k),qim(k)));
                else
                    dzzzm = double(zframe(qjm(k),qim(k))-zframe(qjm(k),qim(k)-1));
                end
            else
                zzzm = 0;
                dzzm = 0;
            end
            
            
            
            % make sure pixel is in sub-region AND within in array bounds
            if    xxxm >= xsubmin && xxxm <= xsubmax ...
                    &&  yyym >= ysubmin && yyym <= ysubmax
                switch pselect
                    case 1
                        rphase=2*pi*double(phaimg(qim(k),qjm(k)));
                    case 5
                        rphase=qqp(k);
                    case {7,9}
                        rphase=qgx(k); % East gradient
                    otherwise
                        error(sprintf('Unknown pselect = %d\n',pselect));
                end
                if abs(rphase) > 0.0 && isfinite(rphase)==1
                    npixinpair = npixinpair+1;
                    ndata1 = ndata1+1;
                    % index of first pixel in this pair
                    if nk == 1
                        ippix1(nk)=1;
                    else
                        %ippix1(nk)=i2+1;
                        ippix1(nk)=ippix2(nk-1)+1;
                    end
                    i1 = ippix1(nk);
                    i2 = ippix1(nk)+npixinpair-1;
                    ippix2(nk) = i2;
                    
                    phao(i2) = rphase;
                    xyzm(1,i2) = xxxm;
                    xyzm(2,i2) = yyym;
                    xyzm(3,i2) = zzzm;
                    dz(i2) = dzzzm;
                    qii1(i2)=qi1(k); % index of quad tree patch
                    qii2(i2)=qi2(k); % index of quad tree patch
                    qiim(i2)=qim(k); % index of quad tree patch midpoint
                    qjj1(i2)=qj1(k); % index of quad tree patch
                    qjj2(i2)=qj2(k); % index of quad tree patch
                    qjjm(i2)=qjm(k); % index of quad tree patch midpoint
                    
                    % number of pixels in patch
                    nnqt(i2) = (qi2(k)-qi1(k)+1) .* (qj2(k)-qj1(k)+1);
                    phasig(i2) = 1./sqrt(nnqt(i2));
                   
                    % Get Latitude, Longitude in degrees
                    if isgeo >= 2 % cartographic
                        %if numel(utmzone0) > 0
                        if isgeo == 2
                            [alat(i2),alon(i2)]=utm2deg(xyzm(1,i2)-x0+xutm0,xyzm(2,i2)-y0+yutm0,utmzone0);
                        elseif isgeo == 3
                            [alat(i2),alon(i2)]=utm2deg(xyzm(1,i2)-x0+xutm0,xyzm(2,i2)-y0+yutm0,utmzone0);
                        else
                            alat(i2) = 0;
                            alon(i2) = 0;
                        end
                    else
%                         alon(i2) = x1 + dc * (c1-1 + double(qjm(k)));
%                         alat(i2) = y1 + dl * (l1-1 + double(qim(k)));
            % 2014-JAN-08 Permute indices
                        alon(i2) = x1 + dc * (c1-1 + double(qim(k)));
                        alat(i2) = y1 + dl * (l1-1 + double(qjm(k)));
                  end
                end
            end
        end % while loop over pixels
        
        if npixinpair < 1
            warning(sprintf('Not enough pixels (%d) in this pair!\n',npixinpair));
        else
            fprintf(1,'Number of good pixels in this pair: %d\n',npixinpair)
            %ndata1 = ndata1 + npixinpair;
            fprintf(1,'Number of observations so far: %d\n',ndata1)
        end
    elseif pselect == 4
        % choose best pixels based on coherence
        clear phase; clear coher; clear elevm; clear xinkm; clear yinkm; clear ikeep; clear jkeep;
        
        coher = cohimg(icenter-halfheight:icenter+halfheight...
            ,jcenter-halfwidth:jcenter+halfwidth); % coherence from 0 to 255
        zframe2 = zframe(icenter-halfheight:icenter+halfheight...
            ,jcenter-halfwidth:jcenter+halfwidth);
        phase2 = phaimg(icenter-halfheight:icenter+halfheight...
            ,jcenter-halfwidth:jcenter+halfwidth);
        phase2 = 2*pi*double(phase2)/256.0; % radians 2010-JAN-11
        disp('Number of pixels in sub-region:');
        npixinsub = numel(coher)
        if numel(zframe2) > npixinsub
            error(sprintf('Problem with dimensions of zframe2\n'));
        end
        if numel(phase2) > npixinsub
            error(sprintf('Problem with dimensions of phase2\n'));
        end
        %  Find the threshold level to get the right number of pixels
        figure
        hist2(reshape(coher,numel(coher),1),50);
        xlabel('Number of occurences');
        ylabel('Coherence');
        disp('Fraction of pixels to retain')
        frac = 1.0-double(npix)/npixinsub
        if frac <= 0
            warning(sprintf('npix (%d) exceeds number of pixels in subregion (%d)\n'...
                ,npix,npixinsub));
            frac = 0.9
        end
        cthresh  = quantile(reshape(coher,numel(coher),1),frac);
        disp('Threshold value for coherence');
        cthresh  = mean(cthresh)
        
        [icoh,jcoh] = find(coher >= cthresh);
        disp('number of good pixels in this pair')
        npixinpair = numel(icoh)
        for k=1:npixinpair
            phase(k) = phase2(icoh(k),jcoh(k)); % wrapped phase in cycles
            elevm(k) = zframe2(icoh(k),jcoh(k)); % DEM elevation in meters
            xinm(k) = xax(jcenter-halfwidth  +jcoh(k));
            yinm(k) = yax(icenter-halfheight +icoh(k));
        end
        disp('number of observations so far')
        ndata1 = ndata1 + npixinpair
        
        % index of first pixel in this pair
        if nk == 1
            ippix1(nk)=1;
        else
            ippix1(nk)=i2+1;
        end
        i1 = ippix1(nk);
        i2 = ippix1(nk)+npixinpair-1;
        
        % values for good pixel
        phao(i1:i2)   = phase;      % wrapped phase in radians
        xyzm(1,i1:i2) = xinm;  % easting in meters
        xyzm(2,i1:i2) = yinm;  % northing in meters
        xyzm(3,i1:i2) = elevm;      % elevation in meters, set to zero!
        dz(i1:i2) = 0; % change in elevation in meters, set to zero!
    elseif ismember(pselect,[0,1,2,3])
        clear phase;
        % 2011-JUN-15 test for zero values
        % and pixels with zero elevation (in water)
        npixinpair = 0;
        for j=1:npix
            %fprintf('\nPixel %5d',j);
            %phase(j) = phaimg(ikeep(j),jkeep(j)); % wrapped phase in DN
            %phase(j) = 2*pi*double(phaimg(ikeep(j),jkeep(j)))/256.0; % wrapped phase in radians
            switch pselect
                case {0,1,2}
                    rphase = 2*pi*double(phaimg(ikeep(j),jkeep(j)))/256.0; % wrapped phase in radians
                case 3
                    rphase = qsp(ikeep(j),jkeep(j)); % wrapped phase in radians, after quad-tree sampling
                otherwise
                    error(sprintf('Unknown pselect = %d\n',pselect));
            end
            if pselect == 0 ...
                    || (abs(rphase) > 0.0 && abs(zframe(ikeep(j),jkeep(j))) > 0)
                kk  = kk+1;
                %fprintf(' kept');
                npixinpair = npixinpair+1;
                phao(kk) = rphase;
                % 2014-JAN-08
                phasig(kk) = 1.0; % no scaling
                % 2012-JAN-12 - finally works!!
                xyzm(1,kk) = xax(jkeep(j));     % easting in meters
                xyzm(2,kk) = yax(ikeep(j));     % northing in meters
                if isgeo < 2
                    alon(kk)   = lonax(jkeep(j)); % longitude in degrees
                    alat(kk)   = latax(jkeep(j)); % latitude in degrees
                end
                %
                %
                %                 if isgeo >= 2 % cartographic
                %                     %   2011-JUN-22
                %                     %                     xyzm(1,kk) = x1 + dc * (c1-1 + double(jkeep(j)));     % easting in meters
                %                     %                     xyzm(2,kk) = y1 + dl * (l1-1 + double(ikeep(j)));     % northing in meters
                %                     %                     % 2012 - JAN - 12
                %                     %                     xyzm(1,kk) = x1 + dc * (double(c1-1) + double(jkeep(j)-1));     % easting in meters
                %                     %                     xyzm(2,kk) = y1 + dl * (double(l1-1) + double(ikeep(j)-1));     % northing in meters
                %                 else % geographic
                %                     %                    xyzm(1,kk) = xstart + dx * double(jkeep(j));     % easting in meters
                %                     %                    xyzm(2,kk) = ystart + dy * double(ikeep(j));     % northing in meters
                %                     %                    % 2011-JUN-22
                %                     %                     xyzm(1,kk) = xstart + dx * (double(jkeep(j)-1));     % easting in meters
                %                     %                     xyzm(2,kk) = ystart + dy * (double(ikeep(j)-1));     % northing in meters
                %                     % %                     % 2012 JAN 11
                %                     %                     xyzm(1,kk) = xstart + dx * (double(jkeep(j)));     % easting in meters
                %                     %                     xyzm(2,kk) = ystart + dy * (double(ikeep(j)));     % northing in meters
                %                     % 2011-JUN-28
                %                     %                     alon(kk) = x1 + dc * (c1-1 + double(jkeep(j))); % longitude in degrees
                %                     %                     alat(kk) = y1 + dl * (l1-1 + double(ikeep(j))); % latitude in degrees
                %                     % 2012-JAN-12
                % %                     xyzm(1,kk) = xstart + dx * (double(jkeep(j)-1));     % easting in meters
                % %                     xyzm(2,kk) = ystart + dy * (double(ikeep(j)-1));     % northing in meters
                % %                     alon(kk)   = x1     + dc * (double(c1-1) + double(jkeep(j)-1)); % longitude in degrees
                % %                     alat(kk)   = y1     + dl * (double(l1-1) + double(ikeep(j)-1)); % latitude in degrees
                %                     xyzm(1,kk) = xax(jkeep(j));     % easting in meters
                %                     xyzm(2,kk) = yax(ikeep(j));     % northing in meters
                %
                %                 end
                xyzm(3,kk) = zframe(ikeep(j),jkeep(j));  % % elevation in meters from DEM
                % topographic relief in east direction
                % change in elevation in meters
                % dz(kk) = 0; % 2011-JUN-22
                if jkeep(j) < mc
                    dz(kk) = double(zframe(ikeep(j),jkeep(j)+1)-zframe(ikeep(j),jkeep(j)));
                else
                    dz(kk) = double(zframe(ikeep(j),jkeep(j))-zframe(ikeep(j),jkeep(j)-1));
                end
            end
        end
        
        disp('Number of data in this pair')
        npixinpair
        
        disp('Number of data so far')
        ndata1 = ndata1 + npixinpair
        
        % index of first pixel in this pair
        if nk == 1
            ippix1(nk)=1;
        else
            ippix1(nk)=i2+1;
        end
        i1 = ippix1(nk);
        i2 = ippix1(nk)+npixinpair-1;
    else
        error(sprintf('Unknown pselect (%d)\n',pselect));
    end % if cases on pselect
    
    %% ************************************************************************
    %% DO THE FOLLOWING TO ALL PIXELS IN PAIR, NO MATTER HOW THEY ARE SELECTED
    %% ************************************************************************
    for ii = i1:i2
        % approximate latitude and longitude from inverse UTM - may have already been done above??
        if isgeo >= 2 && pselect <= 3
            xtemp = xyzm(1,ii)-x0+xutm0;
            ytemp = xyzm(2,ii)-y0+yutm0;
            [alat(ii),alon(ii)]=utm2deg(xtemp,ytemp,utmzone0);
            if ii == i1
                fprintf(1,'xtemp ytemp alat, alon%20.6f %20.6f %20.6f %20.6f %s\n'...
                    ,xtemp,ytemp,alon(ii),alat(ii),utmzone0);
            end
        end
    end
    fprintf(1,'Minimal values in meters         for xyzm %10.1f %10.1f %10.1f \n',min(xyzm(1,i1:i2)),min(xyzm(2,i1:i2)),min(xyzm(3,i1:i2)));
    fprintf(1,'Maximal values in meters         for xyzm %10.1f %10.1f %10.1f \n',max(xyzm(1,i1:i2)),max(xyzm(2,i1:i2)),max(xyzm(3,i1:i2)));
    fprintf(1,'Minimal values in degrees for lon, lat    %10.1f %10.1f        \n',min(alon(i1:i2)),min(alat(i1:i2)));
    fprintf(1,'Maximal values in degrees for lon, lat    %10.1f %10.1f        \n',max(alon(i1:i2)),max(alat(i1:i2)));
    
    % plot image with sample locations
    nf=nf+1;h(nf)=figure;
    imagesc([xsubmin, xsubmax]/1e3,[ysubmax ysubmin]/1e3,double(phaimg(isub,jsub))/256,[-0.5 0.5]);
    hold on;axis xy; axis tight; axis equal;
    %plot(xyzm(1,:)/1000,xyzm(2,:)/1000,'ok');
    plot(xyzm(1,i1:i2)/1000,xyzm(2,i1:i2)/1000,'ok');
    plot(xcenter/1000,ycenter/1000,'+k','MarkerSize',5); % draw a plus sign at center
    plot([xsubmin xsubmax xsubmax xsubmin xsubmin]/1e3...
        ,[ysubmin ysubmin ysubmax ysubmax ysubmin]/1e3...
        ,'k+-');
    xlabel('Easting (km)');ylabel('Northing (km)');
    title(titlestr);
    colorbar;
    cmapblackzero;
    feval(printfun,sprintf('%s_obs_P%02d',runname,i));
    
    %     %plot wrapped phase in histogram
    %     nf=nf+1;h(nf)=figure;hold on;title(titlestr);
    %     %hist([double(phao(i1:i2)) double(phao(i1:i2))+256]/256,64)  %show again, plus 1 cycle to visualize trends
    %     hist2(double(phao(i1:i2))/DNPC,64);
    %     axis xy;%axis ([-0.5 1.5 0 Inf]);
    %     xlabel('wrapped phase (cycles)');ylabel('number of occurences');
    %     feval(printfun,sprintf('%s_Histogram_P%02d',runname,i));
    
    % plot wrapped phase as a function of topographic elevation
    nf=nf+1;h(nf)=figure;hold on;title(titlestr);
    subplot(4,1,1);
    plot(xyzm(1,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k+');
    xlabel('Easting (km)');ylabel('phase change (cycles)');
    subplot(4,1,2);
    plot(xyzm(2,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k+');
    xlabel('Northing (km)');ylabel('phase change (cycles)');
    subplot(4,1,3);
    plot(xyzm(3,i1:i2)/1e3, double(phao(i1:i2))/DNPC,'k+');
    xlabel('topographic elevation (km)');ylabel('phase change (cycles)');
    subplot(4,1,4);
    rdist = sqrt((xyzm(1,i1:i2)-xcenter).^2 +(xyzm(2,i1:i2)-ycenter).^2);
    plot(rdist/1e3, double(phao(i1:i2))/DNPC,'k+');
    xlabel('Radial distance from center of subregion (km)');ylabel('phase change (cycles)');
    feval(printfun,sprintf('%s_PositionVsPhase_P%02d',runname,i));
    
    
    
    if orbits_loaded == 0
        % orbital information varies across scene and thus across pairs
        unitv(1:3,i1:i2) = nan;
        orbvm(1:6,i1:i2) = 0.0; % partial derivative with respect to master orbital parameters
        orbvs(1:6,i1:i2) = 0.0; % partial derivative with respect to slave  orbital parameters
        if numel(orbfile) > 0
            % Calculate orbital quantities at center for master
            [xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum] = readorb(ofnames1{i});
            if min(orbnum) ~= abs(imast(i))
                warning(sprintf('Orbit numbers differ %d %d %d\n',min(orbnum),abs(imast(i)),min(orbnum)-abs(imast(i))));
            end
            [UC, DC, NC, RC, HC, AC, VC, near_mjdC, near_secC] = orbvectors2(loncenter,latcenter,0.0 ...
                ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
            SC = DC/norm(DC); % unit look vector target w.r.t. satellite
            
            for ii=i1:i2
                if ii == i1
                    fprintf(1,'\nCalculating look vectors for MASTER of pair %d for pixels %d to %d from orbit file %s\nPixels numbered:\n',i,i1,i2,ofnames1{i});
                    tstart=tic;
                end
                % calculate unit vector that varies across scene
                if mod(ii,10) == 0; fprintf(1,'%4d ',ii); end
                if mod(ii,100) == 0; fprintf(1,'\n'); end
                %  calculate components of orbit vector
                [U1, D1, N1, R1, H1, A1, V1, near_mjd1, near_sec1] = orbvectors2(alon(ii),alat(ii),xyzm(3,ii)...
                    ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                if abs(near_mjd1 - near_mjdC) > 0
                    fprintf(1,'Longitudes %12.6f %12.6f %12.6f\n',loncenter,alon(ii),loncenter-alon(ii));
                    fprintf(1,'Latitudes  %12.6f %12.6f %12.6f\n',latcenter,alat(ii),latcenter-alat(ii));
                    dt = 86400.0*(near_mjd1 - near_mjdC) + (near_sec1 - near_secC); % time difference in seconds
                    warning(sprintf('different MJD %d %d %d\n',near_mjd1,near_mjdC,near_mjd1 - near_mjdC));
                else
                    dt = near_sec1 - near_secC; % time difference in seconds
                end
                
                unitv(1,ii) = U1(1);    % eastward  component of dimensionless unit vector
                unitv(2,ii) = U1(2);    % northward component of dimensionless unit vector
                unitv(3,ii) = U1(3);    % upward    component of dimensionless unit vector
                
                if orbopt == 1
                    SC = D1/norm(D1); % unit look vector target w.r.t. satellite
                    
                    if ismember(pselect,[7,9]) == 1
                        %  calculate components of baseline
                        [U3, D3, N3, R3, H3, A3, V3, near_mjd3, near_sec3] = orbvectors2(alon(ii)+dlon,alat(ii),xyzm(3,ii)...
                            ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                        S3 = D3/norm(D3); % unit look vector target w.r.t. satellite
                        orbvm(1,ii) = dot(H3,SC) - dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvm(2,ii) = dot(A3,SC) - dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvm(3,ii) = dot(V3,SC)-  dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvm(4,ii) = dt * (dot(H3,SC)-dot(H1,SC)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvm(5,ii) = dt * (dot(A3,SC)-dot(A1,SC)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvm(6,ii) = dt * (dot(V3,SC)-dot(V1,SC)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    else
                        orbvm(1,ii) = dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvm(2,ii) = dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvm(3,ii) = dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvm(4,ii) = dt * dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvm(5,ii) = dt * dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvm(6,ii) = dt * dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    end
                end
            end
            fprintf(1,'\nFinished calculating in %#10.4f seconds\n',toc(tstart));
            
            if orbopt == 1
                % Calculate orbital quantities at center for slave
                [xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum] = readorb(ofnames2{i});
                if min(orbnum) ~= abs(islav(i))
                    warning(sprintf('Orbit numbers differ %d %d %d\n',min(orbnum),abs(islav(i)),min(orbnum)-abs(islav(i))));
                end
                [UC, DC, NC, RC, HC, AC, VC, near_mjdC, near_secC] = orbvectors2(loncenter,latcenter,0.0 ...
                    ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                SC = DC/norm(DC); % unit look vector target w.r.t. satellite
                for ii=i1:i2
                    if ii == i1
                        fprintf(1,'\nCalculating look vectors for SLAVE of pair %d for pixels %d to %d from orbit file %s\nPixels numbered:\n',i,i1,i2,ofnames2{i});
                        tstart=tic;
                    end
                    % calculate unit vector that varies across scene
                    if mod(ii,10) == 0; fprintf(1,'%4d ',ii); end
                    if mod(ii,100) == 0; fprintf(1,'\n'); end
                    %  calculate components of baseline
                    [U2, D2, N2, R2, H2, A2, V2, near_mjd2, near_sec2] = orbvectors2(alon(ii),alat(ii),xyzm(3,ii)...
                        ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                    if abs(near_mjd2 - near_mjdC) > 0
                        fprintf(1,'Longitudes %12.6f %12.6f %12.6f\n',loncenter,alon(ii),loncenter-alon(ii));
                        fprintf(1,'Latitudes  %12.6f %12.6f %12.6f\n',latcenter,alat(ii),latcenter-alat(ii));
                        dt = 86400.0*(near_mjd2 - near_mjdC) + (near_sec2 - near_secC); % time difference in seconds
                        warning(sprintf('different MJD %d %d %d\n',near_mjd2,near_mjdC,near_mjd2 - near_mjdC));
                    else
                        dt = near_sec2 - near_secC; % time difference in seconds
                    end
                    SC = D2/norm(D2); % unit look vector target w.r.t. satellite
                    
                    if ismember(pselect,[7,9]) == 1
                        %  calculate components of baseline
                        [U4, D4, N4, R4, H4, A4, V4, near_mjd4, near_sec4] = orbvectors2(alon(ii)+dlon,alat(ii),xyzm(3,ii)...
                            ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                        S4 = D4/norm(D4); % unit look vector target w.r.t. satellite
                        orbvs(1,ii) = dot(H4,SC) - dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvs(2,ii) = dot(A4,SC) - dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvs(3,ii) = dot(V4,SC)-  dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvs(4,ii) = dt * (dot(H4,SC)-dot(H2,SC)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvs(5,ii) = dt * (dot(A4,SC)-dot(A2,SC)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvs(6,ii) = dt * (dot(V4,SC)-dot(V2,SC)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    else
                        orbvs(1,ii) = dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvs(2,ii) = dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvs(3,ii) = dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvs(4,ii) = dt * dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvs(5,ii) = dt * dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvs(6,ii) = dt * dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    end
                end
            end
            fprintf(1,'\nFinished calculating in %#10.4f seconds\n',toc(tstart));
        else
            fprintf(1,'Assuming constant unit vector [E,N,U] %10.4f %10.4f %10.4f\n'...
                ,unitv0(1),unitv0(2),unitv0(3));
            % constant across scene
            for ii=i1:i2
                unitv(1,ii)=unitv0(1);
                unitv(2,ii)=unitv0(2);
                unitv(3,ii)=unitv0(3);
            end
        end
    end
    
    % close some windows if too many
    if np > 2
        close all
    end
end % loop over files

% Check that orbits are consistent with data
[ndum,nunitv]=size(unitv);
if ndum == 3 && nunitv == ndata1
    orbits_loaded = 1;
    fprintf(1,'Saving orbit information for future use in file named %s\n',orbsavefile);
    if exist('orbvm','var') ~= 1
        orbvm = zeros(6,ndata1);
    end
    if exist('orbvs','var') ~= 1
        orbvs = zeros(6,ndata1);
    end
    save(orbsavefile,'alat','alon','unitv','orbvm','orbvs');
else
    orbits_loaded = 0;
    fprintf(1,'Dimension mismatch number of unit vectors (%d) different from number of pixels(%d)\n',nunitv,ndata1);
    fprintf(1,'Deleting save file %s\n',orbsavefile);
    fdelete(orbsavefile);
    fprintf(1,'Please run again.\n');
    error(sprintf('Dimension mismatch number of unit vectors (%d) different from number of pixels(%d)\n',nunitv,ndata1));
end




% record the names of the files
if nk == np
    if pselect == 6
        fplist = fopen('stackphase.in','w');
        if fplist <= 0
            error('Could not open stackphase.in');
        end
    else
        fplist = [];
    end
    fprintf(1,'Successfully read the following phase files:');
    for i=1:np
        fprintf(1,'%80s\n',pfnames{i});
    end
else
    error(sprintf('Problem(s) reading %d of %d interferogram files\n'...
        ,np-nk,np));
end


disp('Number of pairs with data');np = nk

disp('Number of data');ndata = numel(phao)
disp('Number of data');ndata1

if ndata ~= ndata1
    error(sprintf('NDATA = %d but NDATA1 = %d\n',ndata,ndata1));
end


% i1 = 1;
% i2 = ndata;


% define data types
switch pselect
    case {0,1,2,3,5}
        idatatype = 0;   % phase in radians
    case {7,9}
        idatatype = -1; % eastward phase gradient in radians (per pixel)
    otherwise
        error(sprintf('Undefined pselect = %d\n',pselect));
end

% do not need quad-tree indices
if ismember(pselect,[5,7]) == 0
    qii1 = zeros(ndata,1);
    qii2 = zeros(ndata,1);
    qjj1 = zeros(ndata,1);
    qjj2 = zeros(ndata,1);
end

% 2011-OCT-03 NOTE!!! phao is written in RADIANS to DST file
DST=build_dst(fitfun,xyzm,tepochs,bpest,dops,DD,unitv...
    ,phao,NaN*ones(size(phao)),ippix1,mpercy,idatatype...
    ,dx,dy,dz,orbvm,orbvs,alon,alat...
    ,qii1,qii2,qjj1,qjj2...
    ,phasig);
ierr = write_dst(DST,fitfun,fnamedst);
%check_struct(DST,0); % do not print
check_struct(DST,1); % print min, max, values


% % % plot the pseudoabsolute baselines
% % nf=nf+1;h(nf)=figure;
% % plotbp(tepochs, bpest, DD, species, iuniqorbs, uniqdates, 0);
% % ktour = plotbp(tepochs, bpest, DD, species,iuniqorbs, uniqdates, 3);
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
save;

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));
