function cmdstr = get_executable_name(srcname)
%function cmdstr = get_executable_name(srcname)
% given the name of source code, return the command string required to execute it

% find the proper extension for this computer
[pathstr,fname,ext] = fileparts(srcname);
exeext  = mexext;
exename = strrep(srcname,ext,sprintf('.%s',exeext(4:end)))

%exename = '/Users/feigl/gipht/pha2qls2/pha2qlsV29nowarnings.maci64'
% if numel(getenv('GIPHT_HOME')) > 0
%     cmdstr = strcat(getenv('GIPHT_HOME'),filesep,'src',filesep,exename);
% else
%     cmdstr = strcat('src',filesep,exename);
% end


% if fexist(cmdstr) ~= 1
%     fprintf(1,'Cannot find executable named %s\n',cmdstr);
%     if strfind('.c',ext) > 0
%         fprintf(1,'Consider gcc %s -o %s\n',srcname,exename);
%     end
%     error(sprintf('Cannot find executable named %s\n',cmdstr));
% end

%   exist function returns
%      2 if A is an M-file on MATLAB's search path.  It also returns 2 when
%            A is the full pathname to a file or when A is the name of an
%            ordinary file on MATLAB's search path

if exist(exename,'file') == 2
    cmdstr = which(exename);
else
    fprintf(1,'Cannot find executable named %s\n',exename);
    if strfind('.c',ext) > 0
        fprintf(1,'Consider gcc %s -o %s\n',srcname,exename);
    end
    error(sprintf('Cannot find executable named %s\n',cmdstr));
end

return
