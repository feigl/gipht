function np = stack_pha(pfnames,ncols,nrows,sphnam)
% given an array of names of files containing 1-byte phase files
% make stack of wrapped phase, accumulated over all files
slsnam =sprintf('stack_pha.lst'); % input list of file names

np = numel(pfnames);

% record the names of the files
fplist = fopen(slsnam,'w');
if fplist <= 0
   error(sprintf('Could not open list of phase files %s\n',slsnam));
end
for i=1:np
   fprintf(fplist,'%s\n',pfnames{i});
end
fclose(fplist);

% Make a command line like this:
%/usr1/feigl/GIPhT1.3/src/stack_pha stackpha.in 1520 1540 1 stackpha.pha
% Name of executable, compiled outside of Matlab. DEPENDS ON LOCAL MACHINE !!!
%cmd1='/Users/feigl/Documents/GIPhT1.3/src/stack_pha'; % Kurt's Macbook Pro
%cmd1='/usr1/feigl/GIPhT1.3/src/stack_pha';             % Hengill
cmd1=strrep(which('stack_pha.c'),'.c','');
if fexist(cmd1) <= 0
   fprintf(1,'Cannot find stack_pha executable.\n');
   fprintf(1,'Locate stack_pha.c and compile with:\n');
   fprintf(1,'gcc stack_pha.c -lm -o stack_pha\n');
   error('Cannot find stack_pha executable.');
end

% file names and sizes
cmd2=sprintf('%s %d %d %d %s',slsnam,ncols,nrows,1,sphnam); % 1 byte per pixel
% complete command line
cmd3=sprintf('%s %s',cmd1,cmd2);
fprintf(1,'Starting STACK_PHA with command line:\n%s\n',cmd3);
[unixstat,unixout] = unix(cmd3);
if unixstat == 0
   fprintf(1,'STACK_PHA successful.\n');
else
   error(unixout);
end

% sph=read_pha(sphnam,ncols);
% figure;imagesc(sph);colorbar;cmapblackzero;
% title('wrapped phase stack (256 DN per cycle)')
% xlabel('column index');ylabel('row index');

return

