function VRT = read_vrt(vrt_file_name)
% function VRT = read_vrt(vrt_file_name)
% read VRT header file
% 2021/09/30 Kurt Feigl
%
% INPUT: vrt_file_name with tags:
% <VRTDataset rasterXSize=7200 rasterYSize=3600>
%     <SRS>EPSG:4326</SRS>
%     <GeoTransform>-113.0, 0.0002777777777777778, 0.0, 39.0, 0.0, -0.0002777777777777778</GeoTransform>
%     <VRTRasterBand dataType=Int16 band=1 subClass=VRTRawRasterBand>
%         <SourceFilename relativeToVRT=1>demLat_N38_N39_Lon_W113_W111.dem.wgs84</SourceFilename>
%         <ByteOrder>LSB</ByteOrder>
%         <ImageOffset>0</ImageOffset>
%         <PixelOffset>2</PixelOffset>
%         <LineOffset>14400</LineOffset>
%     </VRTRasterBand>
% </VRTDataset>
% 
% OUTPUTS:
%    VRT is a struct with fields: 
%                 nx: 7200            % number of columns (samples)
%                 ny: 3600            % number of rows (lines)
%                 x0: -113            % X (longitude) coordinate of first pixel
%                 dx: 2.7778e-04      % step size in X coordinate
%                 sx: 0               % scale factor in X
%                 y0: 39              % Y (latitude) coordinate of first pixel
%                 sy: 0               % scale factor in Y
%                 dy: -2.7778e-04     % step size in Y coordinate
%     SourceFilename: 'demLat_N38_N39_Lon_W113_W111.dem.wgs84'


fid=fopen(vrt_file_name,'rt');


while 1
    tline = fgetl(fid);
    if ischar(tline)
        if contains(tline,'<GeoTransform>')
            %disp(tline);
            i1=strfind(tline,'<GeoTransform>')+strlength('<GeoTransform>');
            i2=strfind(tline,'</GeoTransform>')-1;
            %fprintf(1,'%s\n',tline(i1:i2));
            T=sscanf(tline(i1:i2),'%f,%f,%f,%f,%f,%f');
            VRT.x0=T(1);
            VRT.dx=T(2);
            VRT.sx=T(3);
            VRT.y0=T(4);
            VRT.sy=T(5);
            VRT.dy=T(6);
        elseif contains(tline,'<VRTDataset')
            %disp(tline);
            i1=strfind(tline,'rasterXSize=')+strlength('rasterXSize=')+1;
            i2=strfind(tline,'rasterYSize=')-2;
            %fprintf(1,'%s\n',tline(i1:i2));
            VRT.nx=sscanf(tline(i1:i2),'%f');
            i3=strfind(tline,'rasterYSize=')+strlength('rasterYSize=')+1;
            i4=strlength(tline)-2;
            %fprintf(1,'%s\n',tline(i3:i4));
            VRT.ny=sscanf(tline(i3:i4),'%f');
        elseif contains(tline,'<SourceFilename')
            disp(tline);
            i1=strfind(tline,'>')+1;
            i2=strfind(tline,'</SourceFilename>')-1;
            fprintf(1,'%s\n',tline(i1:i2));
            VRT.SourceFilename=sscanf(tline(i1:i2),'%s'); 
        elseif contains(tline,'<VRTRasterBand dataType=')
            %disp(tline);
            i1=strfind(tline,'<VRTRasterBand dataType=')+strlength('<VRTRasterBand dataType=')+1;
            i2=strfind(tline,'" band')-1;
            %fprintf(1,'%s\n',tline(i1:i2));
            VRT.dataType=sscanf(tline(i1:i2),'%s');           
        end
    else
        break
    end   
end
fclose(fid);
return
end

