function ierror = write_oct(name,x)
% WRITE_OCT:	Opens and writes an .oct file (1 unsigned byte per pixel), no header  as from DIAPASON
%
% write_oct('name',x)
%
% Input: 
% name:        A string including the .oct file name
%
% Output: 
% x:           Matrix including the NxM image.
format compact;
fid=fopen(name,'wb');
if fid == -1
  fprintf(1,'Cannot open file %80s\n',name);
  return
else
  fprintf (1,'Opened %s  ',name);
end
count=fwrite(fid,flipud(x)','uint8');
fclose(fid);
if count == numel(x)
   ierror = 0;
else
   ierror = 1;
   error('Count is bad');
end 
return;


