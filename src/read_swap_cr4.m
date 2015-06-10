function z=read_swap_cr4(name,ncols)
%function z=read_swap_cr4(name,ncols)
% Opens and reads an .cr4 file (8  bytes per pixel), no header as from DIAPASON
%
% z=read_cr4('name',ncols)
%
% Input: 
% name:        A string including the .oct file name
% ncols:       The number of columns
%
% Output: 
% x:           Matrix including the NxM image.
% THIS VERSION DOES NOT FLIP
%
%     2009-JAN-16
% Use swapbyte.c for speed
%                  Kurt
%     2013-05-21   generalize for multiple platforms

% 
% %         swap4byte.integer = word4;
% % 
% %         temp0 = swap4byte.character[0];
% %         temp1 = swap4byte.character[1];
% %         swap4byte.character[0] = swap4byte.character[3];
% %         swap4byte.character[1] = swap4byte.character[2];
% %         swap4byte.character[2] = temp1;
% %         swap4byte.character[3] = temp0;
% 
% % read file as a vector
% [c0,count0]=fread(fid0,'char');
% count0
% whos r
% fclose(fid0);
% 
% for i=1:4:count0-3
%     a = c0(i);
%     b = c0(i+1);
%     c = c0(i+2);
%     d = c0(i+3);
%     c1(i) = d;
%     c1(i+1) = c;
%     c1(i+2) = b;
%     c1(i+3) = a;
% end
% 
% fid1=fopen('swap.tmp','w');
% count1=fwrite(fid1,c1,'char');
% count1
% fclose(fid1);
% 
% fid0=fopen(name,'swap.tmp','r');
% read file as a vector

% Swap bytes 4 by 4

% Swap bytes 4 by 4
commandline = sprintf('%s %s 4 swap_cr4.cr4\n',get_executable_name('swapbytes.c'),name);
[status, result] = system(commandline)

if status ~= 0
    error('call to system failed');
end

fid0=fopen('swap_cr4.cr4','r');
if fid0 <= 0
    error(sprintf('Could not open file named %s\n','swap_cr4.cr4'));
end


[r,count]=fread(fid0,'float');
count
%whos r
fclose(fid0);
    
%disp ncols

% split 
N = numel(r);
a = r(1:2:N-1); % real part
b = r(2:2:N);   % imaginary part
z = complex(a,b);
%whos z

% convert vector into an array (aka image)
%mrows = count/ncols;
mrows = numel(z)/ncols;
z=reshape(z,ncols,mrows)';

fprintf (1,'. Read %d rows x %d columns = %ld pixels WITHOUT FLIPPING\n',mrows,ncols,count);
return;


