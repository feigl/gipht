function ierror = write_byt(name,x)
% WRITE_PHA:	Opens and writes an .byt file (1 UNsigned byte per pixel), no header  as from DIAPASON
%
% write_byt('name',x)
%
% Input: 
% name:        A string including the .byt file name
%
% Output: 
% x:           Matrix including the NxM image.
format compact;
fid=fopen(name,'wb');
if fid == -1
  fprintf(1,'Cannot open file %80s\n',name);
  return
% else
%   fprintf (1,'Opened %s  ',name);
end
% 26mar09lap count=fwrite(fid,flipud(x)','int8');
count=fwrite(fid,flipud(rot90(x,1)),'uint8');
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


