function DEM=read_dem(vrt_file_name)
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
VRT=read_vrt(vrt_file_name)
%fid=fopen('demLat_N38_N39_Lon_W113_W111.dem','r')
fid=fopen(VRT.SourceFilename,'r')
%dem=fread(fid,7200*3600,'int16');
DEM=fread(fid,VRT.nx*VRT.ny,lower(VRT.dataType));
size(DEM)
%dem=reshape(dem,3600,7200);
DEM=reshape(DEM,VRT.nx,VRT.ny);
xax=VRT.x0 + VRT.dx*[0:VRT.nx-1];
yax=VRT.y0 + VRT.dy*[0:VRT.ny-1];
if VRT.dy < 0
    fprintf(1,'flipping\n');
    DEM=flipud(DEM);
    yax=fliplr(yax);
end
figure
imagesc(xax,yax,DEM);
axis xy;
axis image;
colormap(jet);
title(VRT.SourceFilename);
colorbar
%printpdf('check_dem.pdf')
printpdf(sprintf(strrep(vrt_file_name,'.vrt','.pdf')));
return
end