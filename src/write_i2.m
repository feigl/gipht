function ierror = write_i2(name,x)
% WRITE_I2:	Opens and writes an .i2 file (2 signed bytes per pixel), no header  as for DIAPASON DEMs
%
% write_i2('name',x)
%
% Input: 
% name:        A string including the .oct file name
%
% Output: 
% x:           Matrix including the NxM image.
% 2012-JAN-03 modified for Matlab R2011b
format compact;
fid=fopen(name,'wb');
if fid == -1
  fprintf(1,'Cannot open file %80s\n',name);
  return
else
  fprintf (1,'Opened %s .',name);
end
count=fwrite(fid,flipud(x)','int16');
fprintf(1,'Wrote %10d elements of 2 bytes each for a total of %10d bytes.\n',count,2*count);

fclose(fid);
%if count == length(x)
if count == numel(x)
   ierror = 0;
else
   ierror = 1;
   error('Count is bad');
end 
return;


