function lut = geotiff_float_to_16bit(infilename,outfilename,map)
%function errcode = geotiff_float_to_16bit(infilename,outfilename)
%
% First, read a geotiff file containing an image array of science values (e.g., radians)
%
% Then, convert the image from science values to an array of indices into a color
% table with 2^16 elements. Such an image is also known as an indexed
% image.
%
% Then, write a geotiff file containing the same geospatial information 
% 20170425 Kurt Feigl

%% initialize
errorcode = 0;

% extract the data from the file
info1 =geotiffinfo(infilename);
[image1,REF1] =geotiffread(infilename);
REF1

%% make a plot if required
debug = true;
if debug
    figure;
    imagesc(info1.BoundingBox(:,1),info1.BoundingBox(:,2),image1);
    xlabel('UTM Easting [m]');
    ylabel('UTM Northing [m]');
    colorbar;
    title(sprintf('Range change in m: %s',infilename),'Interpreter','none');
end

%% scale the values from [min,max] to [0,1]
gray=mat2gray(image1);

%% make a look-up table LUT
%% TODO add NaN values ...
cy=linspace(nanmin(colvec(image1)),nanmax(colvec(image1)),2^16);
ii=1:2^16;
lut = [cy',ii'];

%% convert to indexed image
[indexed, mapgray] = gray2ind(gray,2^16);
whos 
[nrows,ncols] = size(map)
fprintf(1,'first 10 rows of color map\n');
map(1:10,1:3)
fprintf(1,'last 10 rows of color map\n');
map(nrows-10:nrows,1:3)

%% handle missing data
inan = find(isnan(image1) == 1);
indexed(inan) = 1;
nnan = numel(inan);
if  nnan > 0
    fprintf(1,'%s: replacing %d NaN values with index = 1',mfilename,nnan);
end

%% write the Geotiff Image
% geotiffwrite(FILENAME, X, CMAP, R) writes the indexed image in X and
%     its associated colormap, CMAP, to FILENAME. X is spatially referenced
%     by R.

% Force UTM Zone 11
% coordRefSysCode = 32611;
% copy from input file
coordRefSysCode = info1.GeoTIFFCodes.PCS
geotiffwrite(outfilename,indexed,map,REF1,'CoordRefSysCode',coordRefSysCode);

%% verify output
info2 =geotiffinfo(outfilename);
[image2,REF2] =geotiffread(outfilename);
if info2.GeoTIFFCodes.PCS ~= coordRefSysCode
    warning('issue with PCS');
    errorcode = 1;
end

if debug
    figure;
    imagesc(info2.BoundingBox(:,1),info2.BoundingBox(:,2),image2);
    xlabel('UTM Easting [m]');
    ylabel('UTM Northing [m]');
    colorbar;
    title(sprintf('Range change in m: %s',outfilename),'Interpreter','none');
end

return
end
 