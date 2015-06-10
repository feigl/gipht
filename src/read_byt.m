function x=read_byt(name,ncols)
% READ_OCT:	Opens and reads an .oct file (1byte per pixel), no header, as from diapason
%            
%
% x=read_oct('name',ncols)
%
% Input: 
% name:        A string including the .oct file name
% ncols:       The number of columns
%
% Output: 
% x:           Matrix including the NxM image.
format compact;
fid=fopen(name,'r');
if fid == -1
  fprintf(1,'Cannot open file %80s\n',name);
  return
else
  fprintf (1,'Opened %s  ',name);
end
[x,count]=fread(fid,'uint8');
%disp ncols
mrows = count/ncols;
x=reshape(x,ncols,mrows)';
%NO NO NO! x=flipud(x);
fprintf (1,'. Read %d rows x %d columns = %ld bytes WITHOUT FLIPPING\n',mrows,ncols,count);
%disp 'min = ' 
%disp (min(min(x)))
%disp 'max = ' 
%disp (max(max(x)))
return;


