function I = read_aria6(aria6_file_name,mrow,ncol)
%function I = read_aria6(aria6_file_name,mrow,ncol)
% 
% read a file from ARIA
% cat ../S1-IFG_RM_M1S2_TN018_20160117T232659-20170111T232736_s1-poeorb-ec42-v1.2-standard/merged/topophase.flat.geo.hdr
% ENVI
% description = {Data product generated using ISCE}
% samples = 5957
% lines   = 6331
% bands   = 1
% header offset = 0
% file type = ENVI Standard
% data type = 6
% interleave = bip
% byte order = 0
% coordinate system string = {GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137, 298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]]}
% map_info = {Geographic Lat/Lon, 1.0, 1.0, -71.84472222222222, -35.06388888888889, 0.0002777777777777778, 0.0002777777777777778, WGS-84, units=Degrees}

%% BIP format
% https://www.loc.gov/preservation/digital/formats/fdd/fdd000305.shtml
% Identification and description  Explanation of format description terms
% Full name	Band Interleaved by Pixel (BIP) Image Encoding
% Description	
% Band interleaved by pixel (BIP) is one of three primary methods for encoding image data for multiband raster images in the geospatial domain, such as images obtained from satellites. BIP is not in itself an image format, but is a method for encoding the actual pixel values of an image in a file. Images stored in BIP format have the first pixel for all bands in sequential order, followed by the second pixel for all bands, followed by the third pixel for all bands, etc., interleaved up to the number of pixels. The BIP data organization can handle any number of bands, and thus accommodates black and white, grayscale, pseudocolor, true color, and multi-spectral image data.
% 
% Additional information is needed to interpret the image data, such as the numbers of rows, columns, and bands, and relate the image to geospatial locations. This information may be supplied in a file header (typical on the tapes originally used for satellite image data) or in files associated with a raw image data file.
% 
% Relationship to other formats
%     Used by	BIP_file, Band Interleaved by Pixel (BIP) Image File
%     Equivalent to	BIL_enc, Band Interleaved by Line (BIL) image encoding. An alternative ordering for raster image data, a widely used format that offers compromise performance for both spatial and spectral analysis.
%     Equivalent to	BSQ_enc, Band SeQuential (BSQ) image encoding. An alternative ordering for raster image data, optimal for accessing the image spatial information or information for a particular spectral band.

% http://earthdef.caltech.edu/boards/4/topics/578?r=579#message-579
%  RE: Interferogram to deformation map - Added by Eric Fielding 11 months ago
% 
% Hi Katleen,
% 
% I see now that you are trying to display the wrapped interferogram (.flat.geo), not the unwrapped interferogram (.unw.geo). The wrapped interferogram is complex values, which the MDX program knows how to display but most GIS programs do not understand. The MDX program has extra code that takes the complex values "Type=CFloat32" (two bands of real and imaginary components band-interleaved by pixel) and converts them to phase and magnitude.
% 
% To work with the interferograms in GIS programs you should add the unwrapping step to your processing. I assume you are using insarApp. You can find these lines in the example/input_files/insarApp.xml:
% 
%        <!--<property name="unwrap">True</property>-->
%         <!--<property name="unwrapper name">snaphu_mcf</property>-->
% 
% You should uncomment the lines and put them in your insarApp.xml file within the main block of properties.
% 
% ++Eric

%% open and read the file
%fid=fopen(aria6_file_name,'r','ieee-be'); % big endian
fid=fopen(aria6_file_name,'r','ieee-le'); % little endian
%[I,count1]=fread(fid,[ncol,mrow],'float32');
[I,count1]=fread(fid,[ncol,mrow],'int32');

fprintf(1,'Number of 4-byte numbers read     = %ld\n',count1);

fprintf(1,'Number of pixels         expected = %ld\n',mrow * ncol); 
fprintf(1,'Number of pixels             read = %ld\n',count1);

return

end

