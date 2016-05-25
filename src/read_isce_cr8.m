function [a,b] = read_isce_cr8(file_name,nrows,ncols)
%function [amp,pha] = read_bil(bil_file_name,mrow,ncol)
% Read an interferogram from ISCE in interleaved format
% 4 bytes for real part
% 4 bytes for imaginary part
% -rw-rw-r--+ 1 feigl  15  10579968 Nov 24 12:10 filt_topophase.flat.geo
% -rw-rw-r--+ 1 feigl  15      3091 Nov 24 14:38 filt_topophase.flat.geo.xml

% Kurt Feigl & Elena Baluyut 20151124 


%% open the file
fid=fopen(file_name,'r','ieee-le'); % little endian

%% read the file in column order
[r1,count1]=fread(fid,'real*4');

fprintf(1,'Number of 4-byte numbers read     = %ld\n',count1);

fprintf(1,'Number of pixels         expected = %ld\n',nrows * ncols); 
fprintf(1,'Number of pixels             read = %ld\n',count1);

r1 = reshape(r1,ncols*2,nrows);
r1 = transpose(r1);


%% parse
a = r1(1:nrows,1:2:ncols*2); % real part
b = r1(1:nrows,2:2:ncols*2); % imaginary part

end

