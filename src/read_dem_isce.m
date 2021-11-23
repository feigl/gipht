function DEM=read_dem_isce(vrt_file_name)
% read and verify DEM in ISCE format, as specified by VRT vrt_file_name
% 2021/09/30 Kurt Feigl

% verify dem
% <VRTDataset rasterXSize="7200" rasterYSize="3600">
%     <SRS>EPSG:4326</SRS>
%     <GeoTransform>-113.0, 0.0002777777777777778, 0.0, 39.0, 0.0, -0.0002777777777777778</GeoTransform>
%     <VRTRasterBand dataType="Int16" band="1" subClass="VRTRawRasterBand">
%         <SourceFilename relativeToVRT="1">demLat_N38_N39_Lon_W113_W111.dem</SourceFilename>
%         <ByteOrder>LSB</ByteOrder>
%         <ImageOffset>0</ImageOffset>
%         <PixelOffset>2</PixelOffset>
%         <LineOffset>14400</LineOffset>
%     </VRTRasterBand>
% </VRTDataset>

%VRT=read_vrt('demLat_N38_N39_Lon_W113_W111.dem.wgs84.vrt')
DEM=read_vrt(vrt_file_name)
%fid=fopen('demLat_N38_N39_Lon_W113_W111.dem','r')
fid=fopen(DEM.SourceFilename,'r');
%dem=fread(fid,7200*3600,'int16');
DEM.h=fread(fid,DEM.nx*DEM.ny,lower(DEM.dataType));
% check for zero or negative values
izero = find(DEM.h <= 0);
if numel(izero) > 0
    DEM.h(izero)=NaN;
    warning(sprintf('Found %d values less than or equal to zero. Setting to NaN\n'));
end   

figure;
histogram(colvec(DEM.h));
xlabel('DEM value');
ylabel('count');
title(DEM.SourceFilename,'Interpreter','none');
printpdf(sprintf(strrep(vrt_file_name,'.vrt','.histogram.pdf')));


size(DEM.h)
%% reshape row-wise
% A=[1,2,3,4,5,6]'
% 
% A =
% 
%      1
%      2
%      3
%      4
%      5
%      6

%transpose(reshape(A,3,2))
% ans =
% 
%      1     2     3
%      4     5     6
DEM.h=transpose(reshape(DEM.h,DEM.nx,DEM.ny));

%% construct vector axes
DEM.xax=DEM.x0 + DEM.dx*[0:DEM.nx-1];  % row vector
DEM.yax=DEM.y0 + DEM.dy*[0:DEM.ny-1]'; % column vector

%% flip if necessary
if DEM.dy < 0
    fprintf(1,'flipping\n');
    DEM.h=flipud(DEM.h);
    DEM.yax=flipud(DEM.yax);
end
%% flop if necessary
if DEM.dx < 0
    fprintf(1,'flopping\n');
    DEM.h=fliplr(DEM.h);
    DEM.xax=fliplr(DEM.xax);
end
figure
imagesc(DEM.xax,DEM.yax,DEM.h);
colormap(jet);
axis xy;
axis image;
xlabel('Longitude [deg E]');
ylabel('Latitude [deg N]');
C =colorbar;
C.Label.String = 'Elevation [m]';
title(DEM.SourceFilename,'Interpreter','none');
printpdf(sprintf(strrep(vrt_file_name,'.vrt','.map.pdf')));

return
end