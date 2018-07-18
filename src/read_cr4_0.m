function z = read_cr4(cr4_file_name,mrow,ncol)
% Read a file in complex real*4 format
% 20180709 Kurt Feigl


%% open and read the file
fid=fopen(cr4_file_name,'r','ieee-le'); % little endian
% Helene
[r1,count1]=fread(fid,[ncol,mrow*2],'real*4');
% Kurt
%[r1,count1]=fread(fid,[mrow,Inf],'float32');

fprintf(1,'Number of 4-byte numbers read     = %ld\n',count1);

fprintf(1,'Number of pixels         expected = %ld\n',mrow * ncol); 
fprintf(1,'Number of pixels             read = %ld\n',count1/4/2);

%% parse
%pha = (r(1:ncols,2:2:mrows*2))'; % phase
% Helene
% amp=(r1(1:ncol,1:2:mrow*2))'; % amplitude
% pha = (r1(1:ncol,2:2:mrow*2))'; % phase
% Kurt
% amp=(r1(1:mrow,1:2:ncol*2))'; % amplitude
% pha = (r1(1:mrow,2:2:ncol*2))'; % phase
% get complex
% x = (r1(1:ncol,1:2:mrow*2))'; % real part
% y = (r1(1:ncol,2:2:mrow*2))'; % imaginary part
x = (r1(1:ncol,1:2:mrow*2))'; % real part
y = (r1(1:ncol,2:2:mrow*2))'; % imaginary part
z = complex(x,y);
return;


end

