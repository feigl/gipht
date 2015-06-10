function DST = read_dst(fnamedat)
%function DST = read_dst(fnamedat)
% read a DST structure from a flat ASCII file named fnamedat

% open file for data
fid = fopen(fnamedat,'r');
if fid <= 0
    error(sprintf('Cannot open data file called %s\n',fnamedat));
end

% number of columns
%ncols = 19;
[ncols, nread] = fscanf(fid,'%d\n',1);
fprintf(1,'Number of columns in file is %d\n',ncols);

% count the number of lines
nlines = -2;
aline = ' ';
while aline ~= -1 
    %[A, nread] = fscanf(fid,fmt,ncols);
    aline = fgetl(fid);
    if aline ~= -1
       nlines = nlines+1;
    end
end
fprintf(1,'Number of pixels in file is %d\n',nlines);
frewind(fid);


% write a format statement
fmt = '%g'; % 
% must match number of real variables in fprintf statement below
for i=1:ncols-1
    fmt = sprintf('%s %s',fmt,'%g');
end
fmt = sprintf('%s\\n',fmt);

% initialize structure for speed
DST.i               = zeros(nlines,1); % index
DST.idatatype       = zeros(nlines,1); % type of data
DST.k               = zeros(nlines,1); % index to pair
DST.kmast           = zeros(nlines,1); % index to master
DST.kslav           = zeros(nlines,1); % index to slave
DST.phaobs          = zeros(nlines,1); % observed phase in radians
DST.phamod          = zeros(nlines,1); % modeled phase in radians
DST.tmast           = zeros(nlines,1); % master time (epoch) in decimal years
DST.tslav           = zeros(nlines,1); %  slave time (epoch) in decimal years
DST.x               = zeros(nlines,1); % easting coordinate of obs point in meters
DST.y               = zeros(nlines,1); % northing coordinate of obs point in meters
DST.z               = zeros(nlines,1); % upping coordinate of obs point in meters
DST.uvx             = zeros(nlines,1); % east  component of unit vector pointing from ground to sat
DST.uvy             = zeros(nlines,1); % north component of unit vector pointing from ground to sat
DST.uvz             = zeros(nlines,1); % up    component of unit vector pointing from ground to sat
DST.mpercy          = zeros(nlines,1); % range change per fringe in meters per cycle
DST.bmast           = zeros(nlines,1); % BPerp for master (meters)
DST.bslav           = zeros(nlines,1); % BPerp for slave  (meters)
DST.dmast           = zeros(nlines,1); % Doppler for master (PRF or Hz?)
DST.dslav           = zeros(nlines,1); % Doppler for slave (PRF or Hz?)
DST.x0              = zeros(nlines,1); % reference  easting coordinate in meters
DST.y0              = zeros(nlines,1); % reference northing coordinate in meters
DST.z0              = zeros(nlines,1); % reference vertical coordinate in meters
DST.dx              = zeros(nlines,1); % pixel or patch dimension in  easting direction in meters
DST.dy              = zeros(nlines,1); % pixel or patch dimension in northing direction in meters
DST.dz              = zeros(nlines,1); % vertical step in eastward direction in meters
DST.qii1            = zeros(nlines,1); % quadtree indices
DST.qii2            = zeros(nlines,1); % 
DST.qjj1            = zeros(nlines,1); % 
DST.qjj2            = zeros(nlines,1); % 
DST.phasig          = zeros(nlines,1); % scale factor for measurement uncertainty


k = 0;
nread = ncols;
aline = fgetl(fid); % skip header line
aline = fgetl(fid); % skip header line
aline = ' ';
while aline ~= -1
    %[A, nread] = fscanf(fid,fmt,ncols);
    aline = fgetl(fid);
    if aline ~= -1
        [A, nread] = sscanf(aline, fmt, ncols);
        if nread == ncols
            tstart = tic;
            k = k+1;
            j = 1;
            DST.i(k)            = A(j);j=j+1; % index to datum
            DST.idatatype(k)    = A(j);j=j+1; % type of data
            DST.k(k)            = A(j);j=j+1; % index to pair
            DST.kmast(k)        = A(j);j=j+1; % index to master
            DST.kslav(k)        = A(j);j=j+1; % index to slave
            DST.phaobs(k)       = A(j);j=j+1; % observed phase in radians
            DST.phamod(k)       = A(j);j=j+1; % modeled phase in radians
            DST.tmast(k)        = A(j);j=j+1; % master time (epoch) in decimal years
            DST.tslav(k)        = A(j);j=j+1; %  slave time (epoch) in decimal years
            DST.x(k)            = A(j);j=j+1; % easting coordinate of obs point in meters
            DST.y(k)            = A(j);j=j+1; % northing coordinate of obs point in meters
            DST.z(k)            = A(j);j=j+1; % upping coordinate of obs point in meters
            DST.uvx(k)          = A(j);j=j+1; % east  component of unit vector pointing from ground to sat
            DST.uvy(k)          = A(j);j=j+1; % north component of unit vector pointing from ground to sat
            DST.uvz(k)          = A(j);j=j+1; % up    component of unit vector pointing from ground to sat
            DST.mpercy(k)       = A(j);j=j+1; % range change per fringe in meters per cycle
            DST.bmast(k)        = A(j);j=j+1; % BPerp for master (meters)
            DST.bslav(k)        = A(j);j=j+1; % BPerp for slave  (meters)
            DST.dmast(k)        = A(j);j=j+1; % Doppler for master (PRF or Hz?)
            DST.dslav(k)        = A(j);j=j+1; % Doppler for slave (PRF or Hz?)
            DST.x0(k)           = A(j);j=j+1; % reference  easting coordinate in meters
            DST.y0(k)           = A(j);j=j+1; % reference northing coordinate in meters
            DST.z0(k)           = A(j);j=j+1; % reference vertical coordinate in meters
            DST.dx(k)           = A(j);j=j+1; % pixel or patch dimension in  easting direction in meters
            DST.dy(k)           = A(j);j=j+1; % pixel or patch dimension in northing direction in meters
            DST.dz(k)           = A(j);j=j+1; % vert relief in meters
            DST.alond(k) = A(j);j=j+1;
            DST.alatd(k) = A(j);j=j+1;
            DST.orbm1(k) = A(j);j=j+1;
            DST.orbm2(k) = A(j);j=j+1;
            DST.orbm3(k) = A(j);j=j+1;
            DST.orbm4(k) = A(j);j=j+1;
            DST.orbm5(k) = A(j);j=j+1;
            DST.orbm6(k) = A(j);j=j+1;
            DST.orbs1(k) = A(j);j=j+1;
            DST.orbs2(k) = A(j);j=j+1;
            DST.orbs3(k) = A(j);j=j+1;
            DST.orbs4(k) = A(j);j=j+1;
            DST.orbs5(k) = A(j);j=j+1;
            DST.orbs6(k) = A(j);j=j+1;
            
            DST.qii1     = A(j);j=j+1;
            DST.qii2     = A(j);j=j+1;
            DST.qjj1     = A(j);j=j+1;
            DST.qjj2     = A(j);j=j+1;
            DST.phasig   = A(j);j=j+1;
        end
        if j-1 ~= ncols
            error(sprintf('Miscount. j-1 (%d) does not equal ncols (%d)\n',j-1,ncols));
        end
        if mod(k,10000) == 0
            fprintf(1,'Read line %5d %12.4g\n',k,toc(tstart));
        end
    end
end
fclose(fid);


%
%
%             ,i,k,kmast,kslav...
%             ,xd(i),phamod,tepochs(kmast),tepochs(kslav)...
%             ,xyzm(1,i),xyzm(2,i),xyzm(3,i)...
%             ,unitv(1,i),unitv(2,i),unitv(3,i)...
%             ,mpercy...
%             ,bpest(kmast),bpest(kslav)...
%             ,dops(kmast),dops(kslav));
%
% % loop file
% kount = 0;
% for k=1:np
%     i1 = ippix1(k);
%     if k < np
%         %i2 = ippix1(k+1);
%         i2 = ippix1(k+1)-1;
%     else
%         i2 = ndata;
%     end
%     npixinpair=i2-i1+1;
%
%     kmast = find(DD(k,:) == -1);
%     kslav = find(DD(k,:) == +1);
%
%     for i=i1:i2
%         kount = kount+1;
%         fprintf(fid,fmt...
%             ,i,k,kmast,kslav...
%             ,xd(i),phamod,tepochs(kmast),tepochs(kslav)...
%             ,xyzm(1,i),xyzm(2,i),xyzm(3,i)...
%             ,unitv(1,i),unitv(2,i),unitv(3,i)...
%             ,mpercy...
%             ,bpest(kmast),bpest(kslav)...
%             ,dops(kmast),dops(kslav));
%     end
% end
%
% % check count
% if kount == ndata
%     fprintf(1,'Wrote meta data for %d pixels to file %s\n',kount,fnamedat);
%     fclose(fid);
%     ierr = 0;
% else
%     error(sprintf('kount (%d) not equal to ndata (%d)\n',kount,ndata));
%     ierr = 1;
% end

return

