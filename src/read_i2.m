function x=read_i2(name,ncols)

% READ_I2:	Opens and reads an .i2 file (2 bytes per pixel), no header
%               as from DIAPASON
%
% x=read_i2('name',ncols)
%
% Input: 
% name:        A string including the .i2 file name
% ncols:       The number of columns
%
% Output: 
% x:           Matrix including the NxM image.
% THIS VERSION DOES NOT FLIP
%
%     2009-MAY-10
%                  Kurt

fid=fopen(name,'r');
if fid == -1
  error(sprintf('Cannot open file %80s\n',name));
else
  fprintf (1,'Opened %s',name);
end


%[x,count]=fread(fid,'int16');
[x,count]=fread(fid,'int16=>int16');
%[x,count]=fread(fid,'short');
%[x,count]=fread(fid,'uint16');

%count
%ncols
mrows = count/ncols;
 
x=reshape(x,ncols,mrows)'; % This is a double operation, but correct. Kurt 2009-APR-01
%x=reshape(x,mrows,ncols);  % This is NOT equivalent!! It is wrong.

%%%% NO NO NO x=flipud(x);
%disp 'min = ' 
%disp (min(min(x)))
%disp 'max = ' 
%disp (max(max(x)))
fprintf (1,'. Read 2 x %d rows x %d columns = %ld bytes WITHOUT FLIPPING as int16 \n',mrows,ncols,2*count);
return;


