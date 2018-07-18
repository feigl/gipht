function z = read_cr4(cr4_file_name,nrows,mcols)
% Read a file in complex real*4 format
% 20180709 Kurt Feigl

% https://www.mathworks.com/help/matlab/ref/fread.html
% Read all the data in the file into a vector of class double. By default, fread reads a file 1 byte at a time, interprets each byte as an 8-bit unsigned integer (uint8), and returns a double array.
% 
% fileID = fopen('nine.bin');
% A = fread(fileID)
% 
% A = 9×1
% 
%      1
%      2
%      3
%      4
%      5
%      6
%      7
%      8
%      9
% 
% fread returns a column vector, with one element for each byte in the file. 
% Read the first six values into a 3-by-2 array. Specify that the source data is class uint16.
% 
% fileID = fopen('nine.bin');
% A = fread(fileID,[3,2],'uint16')
% 
% A = 3×2
% 
%      1     4
%      2     5
%      3     6
% 
% fread returns an array populated column-wise with the first six values from the file, nine.bin.
% 
% Return to the beginning of the file.
% 
% frewind(fileID)
% 
% Read two values at a time, and skip one value before reading the next values. Specify this format using the precision value, '2*uint16'. Because the data is class uint16, one value is represented by 2 bytes. Therefore, specify the skip argument as 2.
% 
% precision = '2*uint16';
% skip = 2;
% B = fread(fileID,[2,3],precision,skip)
% 
% B = 2×3
% 
%      1     4     7
%      2     5     8
% 
% fread returns a 2-by-3 array populated column-wise with the values from nine.bin. 

%% open and read the file
fid=fopen(cr4_file_name,'r','ieee-le'); % little endian
% 
nrows
mcols
npixels = 2*nrows*mcols
[r,count1]=fread(fid,npixels,'float32');

fprintf(1,'Number of 4-byte numbers read     = %ld\n',count1);

fprintf(1,'Number of pixels         expected = %ld\n',nrows * mcols); 
fprintf(1,'Number of pixels             read = %ld\n',count1/2);

%% parse
x = r(1:2:npixels-1); % real part
y = r(2:2:npixels); % imaginary part

z = complex(x,y);
z = reshape(z,mcols,nrows);
z = z';
return;


end

