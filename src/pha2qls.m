function npatches = pha2qls(sphnam,ncols,nrows,qphnam,grxnam,grynam,qlsnam,ithresh,minpix,maxcmd,maxpix,pha2qlsname)  
% perform quadtree partitioning using pha2qls3 program 
% that is sensitive to slope
% 2011-07-04

npatches = NaN;
% Make a command line like this:
%/usr1/feigl/GIPhT1.3/src/pha2qls x_pstackpha.pha 1520 1540 x_qstackpha.pha x_qpha2qls.i2 32 9
% Name of executable, compiled outside of Matlab. DEPENDS ON LOCAL MACHINE !!!
%cmd1='/Users/feigl/Documents/GIPhT1.3/src/pha2qls'; % Kurt's Macbook Pro
%cmd1='/usr1/feigl/GIPhT1.3/src/pha2qls';             % Hengill
%cmd1=strrep(which('pha2qls3.c'),'.c','.exe');
% ../src/pha2qls psp_test.pha 501 501 -P qsp_test.pha -G -L 16 -N 4 -M 32

% srcname = 'pha2qls.c';
% exeext  = mexext;
% exename = strrep(srcname,'.c',sprintf('.%s',exeext(4:end)));
% 
% % if exist('which') == 5
% %     % NOTE: COMMAND "WHICH" DOES NOT APPEAR TO WORK IN COMPILED VERSION
% %     fprintf(1,'Using "which" command to search for for executable named %s\n',exename);
% %     cmd1 = which(exename)
% %     if numel(cmd1) > 0
% %         fprintf(1,'Found executable named %s %s\n',exename,cmd1);
% %     else
% %         fprintf(1,'Command "which" failed while looking for executable named %s\n',exename);
% %         error('Which failed');
% %     end
% % else
% %     fprintf(1,'Cannot use "which" command to look for executable named %s\n',exename);
% %     cmd1 = exename;
% %     warning('Cannot find which command');
% % end
% %if fexist(exename) == 1
% if fexist(strcat('.',filesep,exename)) == 1
%     cmd1 = strcat('.',filespe,exename)
% % elseif fexist(pha2qlsname) == 1
% %     %cmd1 = pha2qlsname
% %     cmd1 = strcat(pha2qlsname,sprintf('.%s',exeext(4:end)))
% else
%     if numel(getenv('GIPHT_HOME')) > 0
%         cmd1 = strcat(getenv('GIPHT_HOME'),filesep,'src',filesep,exename)
%     end
% end
% 
% if fexist(cmd1) ~= 1
%     fprintf(1,'Cannot find executable named %s\n',cmd1);
%     error(sprintf('Cannot find executable'));
% end

cmd1 = get_executable_name('pha2qls.c');

% file names and sizes, ithresh = 4 DN, minpix = 16
cmd2=sprintf('%s %d %d -V -P %s -X %s -Y %s -L %d -N %d -M %d -Q %d'...
    ,sphnam,ncols,nrows,qphnam...
    ,grxnam,grynam ...
    ,ithresh,minpix,maxcmd,maxpix);
% complete command line
cmd3=sprintf('%s %s',cmd1,cmd2);
fprintf(1,'Starting pha2qls with command line:\n%s\n',cmd3);
%[unixstat,unixout] = unix(cmd3);
[unixstat,unixout] = system(cmd3);
if unixstat == 0
   fprintf(1,'pha2qls successful.\n');
   %whos unixout
   %unixout
   key = 'N(OK patches)  =';
   k = strfind(unixout,key)+numel(key)+1;
   if k > 0
      npatches = str2num(unixout(k:k+6));
   end
else
   error(['Error executing PHA2QLS ' unixout]);
end
return

   
   
   
