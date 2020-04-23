

DR.filename = 'drhomaskd_utm.tif'
DR.info  =geotiffinfo(DR.filename);
[DR.image,DR.R] =geotiffread(DR.filename)
imagesc(DR.info.BoundingBox(:,1),DR.info.BoundingBox(:,2),DR.image);
xlabel('UTM Easting [m]');
ylabel('UTM Northing [m]');
colorbar;
title(sprintf('Range change in m: %s',DR.filename),'Interpreter','none');

% scale the values
DR.gray=mat2gray(DR.image);

% convert to index
[DR.indexed, DR.map] = gray2ind(DR.gray,2^16);


% geotiffwrite(FILENAME, X, CMAP, R) writes the indexed image in X and
%     its associated colormap, CMAP, to FILENAME. X is spatially referenced
%     by R.

outfilename = strrep(DR.filename,'.tif','_ind.tif')
% geotiffwrite(outfilename,DR.indexed,DR.map,DR.R,'GeoKeyDirectoryTag',DR.info.GeoTIFFTags)

%   grdconvert $IN.grd=nf tmpout.tif=gd:GTiFF
%    # zone 11
%    #gdal_translate -a_srs EPSG:32611 tmpout.tif $IN.tif
%    # any zone
%    gdal_translate -a_srs EPSG:326$UTMZONE tmpout.tif $IN.tif
coordRefSysCode = 32611;
geotiffwrite(outfilename,DR.indexed,DR.map,DR.R,'CoordRefSysCode',coordRefSysCode);
 