function ierror = write_pha(name,x)
% WRITE_PHA:	Opens and writes an .pha file (1 signed byte per pixel), no header  as from DIAPASON
%
% write_pha('name',x)
%
% Input: 
% name:        A string including the .pha file name
%
% Output: 
% x:           Matrix including the NxM image.
% 2011-JUN-21 translate NaN to null
format compact;
inan = find(isfinite(x)==0);
x(inan)=0;
fid=fopen(name,'wb');
if fid == -1
  error(sprintf('Cannot open file %80s\n',name));
end
% 26mar09lap count=fwrite(fid,flipud(x)','int8');
count=fwrite(fid,flipud(rot90(x,1)),'int8');
[m,n] = size(x);
fprintf(1,'Wrote %s with %8d pixels of 1 byte each with %d lines by %d columns.\n',name,count,m,n);
fclose(fid);
if count == m*n
   ierror = 0;
else
   ierror = 1;
   error('Count is bad');
end 
return;


