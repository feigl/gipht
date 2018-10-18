function [b1,b2] = read_bil2r4(bil2r4_file_name,nrows,mcols)
% Read a BIL file with two channels, each real*4
% 20180912 Kurt Feigl

nrows;
mcols;
nbands = 2;

npixels = nrows * mcols;
nwords = nbands*npixels;



%% Binary Interleaved by Line
% http://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/bil-bip-and-bsq-raster-files.htm
% Band interleaved by line data stores pixel information band by band for
% each line, or row, of the image. For example, given a three-band image,
% all three bands of data are written for row 1, all three bands of data
% are written for row 2, and so on, until the total number of rows in the
% image is reached. The following diagram illustrates BIL data for a
% three-band dataset

% Three-element vector of integers consisting of [height, width, N], where
% height is the total number of rows 
% width is the total number of elements in each row 
% N is the total number of bands.

% This will be the dimensions of the data if it is read in its entirety.

% Read the data using Band-Interleaved-by-Line format.
% 
% im3 = multibandread(filename, [rows cols bands], ...
%                     'double', 0, 'bil', 'ieee-le')

r = multibandread(bil2r4_file_name, [nrows mcols nbands], ...
                     'float32', 0, 'bil', 'ieee-le'); 
b1 = r(:,:,1);
b2 = r(:,:,2);

return;

end

