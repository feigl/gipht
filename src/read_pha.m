function pha=read_pha(name,ncols,nrows)
% READ_PHA:	Opens and reads a .pha file (1byte per pixel), no header, as from DIAPASON
%
% pha=read_pha('name',ncols)
%
% Input: 
% name:        A string including the .oct file name
% ncols:       The number of columns
%
% Output: 
% pha:           Matrix including the NxM image, values [-128,127]
%
% 2009-APR-01 Kurt THIS VERSION DOES NOT FLIP
% 2009-JUN-18 Kurt return int8
% 2011-JUL-20 pad or truncate if necessary
%     
%  
format compact;
fid=fopen(name,'r');
if fid == -1
  fprintf(1,'Cannot open file %80s\n',name);
  return
else
  fprintf (1,'Opened %25s\n',name);
end
%[pha,count]=fread(fid,'int8');
[pha,count]=fread(fid,'int8=>int8');
%disp ncols
if exist('nrows','var') == 0
   nrows = count/ncols;
end
if abs(mod(nrows,1)) > 0 
    warning(sprintf('Number of rows (%12.2f) is not integer. Padding with ceil(). \n',nrows));
    nrows = ceil(nrows);
end

if count ~= nrows*ncols
    warning(sprintf('Number of bytes (%d) in phase file %s is not equal to product of nrows (%d) and ncols (%d) = (%d)'...
        ,count,name,nrows,ncols,nrows*ncols));
    if count < nrows*ncols
        warning(sprintf('Paddding %d pixels with zeros\n',nrows*ncols-count));
        pha(count+1:nrows*ncols) = 0;
    else
        warning(sprintf('Truncating %d pixels\n',count-nrows*ncols));
        pha = pha(1:nrows*ncols);
    end
end

fprintf (1,'Reshaping to %d rows pha %d columns = %ld bytes as signed 1-byte integers\n',nrows,ncols,count);

pha=reshape(pha,ncols,nrows)'; % This is a double operation, but correct. Kurt 2009-APR-01
%pha=reshape(pha,nrows,ncols);  % This is NOT equivalent!! It is wrong.

%%NO NO NO pha=flipud(pha);
%fprintf (1,'. Read %d rows pha %d columns = %ld bytes WITHOUT FLIPPING\n',nrows,ncols,count);
fprintf (1,'Read %d rows pha %d columns = %ld bytes as signed 1-byte integers\n',nrows,ncols,count);
% disp 'min = ' 
% disp (min(min(pha)))
% disp 'max = ' 
% disp (max(max(pha)))
return;


