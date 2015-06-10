function x=read_swap_r4(name,ncols)
%function x=read_swap_r4(name,ncols)
% Opens and reads an .r4 file (2 bytes per pixel), no header as from DIAPASON
%
% Input: 
% name:        A string including the .oct file name
% ncols:       The number of columns
%
% Output: 
% x:           Matrix including the NxM image.
% THIS VERSION DOES NOT FLIP
%
%     2009-JAN-16  Kurt
%     2013-05-21   generalize for multiple platforms

fid=fopen(name,'r');
if fid == -1
  error(sprintf(1,'Cannot open file %80s\n',name));
else
  fprintf (1,'Opened %s  ',name);
  fclose(fid);
end

% Swap bytes 4 by 4
commandline = sprintf('%s %s 4 swap_r4.r4\n',get_executable_name('swapbytes.c'),name);

[status, result] = system(commandline)

% This is painfully slow
% % read file as a vector
% [c0,count0]=fread(fid,'char');
% count0
% whos r
% fclose(fid);
% 
% fprintf(1,'Swapping bytes by fours, slowly...\n');
% for i=1:4:count0-3
%     a = c0(i);
%     b = c0(i+1);
%     c = c0(i+2);
%     d = c0(i+3);
%     c1(i) = d;
%     c1(i+1) = c;
%     c1(i+2) = b;
%     c1(i+3) = a;
%     progress = 100 * i/count0;
%     if mod(i,1e4) == 1
%         fprintf (1,'%.3f percent done\n',progress);
%     end
% end
% fprintf(1,'\nDone\n');
% 
% fid=fopen('swap.tmp','w');
% count1=fwrite(fid,c1,'char');
% count1
% fclose(fid);

fid=fopen('swap_r4.r4','r');
[x,count]=fread(fid,'float');
%disp ncols
mrows = count/ncols;
x=reshape(x,ncols,mrows)';
%%% NO NO NO x=flipud(x);
%disp 'min = ' 
%disp (min(min(x)))
%disp 'max = ' 
%disp (max(max(x)))
fprintf (1,'. Read %d rows x %d columns = %ld bytes WITHOUT FLIPPING\n',mrows,ncols,count);
return;


