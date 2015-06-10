function giphtpath
% set path for GIPHT
home = getenv('HOME');
giphtv = 'gipht';
gipht_home = strcat(home,filesep,giphtv)
if exist(gipht_home,'dir') == 7
    fprintf(1,'Setting environment variable GIPHT_HOME  to %s\n',gipht_home);
    setenv('GIPHT_HOME',gipht_home);
    fprintf(1,'Environment variable GIPHT_HOME is now set to %s\n',getenv('GIPHT_HOME'));
    p = strcat(gipht_home,filesep,'src',pathsep,gipht_home,filesep,'extern',pathsep);
    addpath(p,'-BEGIN');
    fprintf(1,'Matlab command search path is now: %s\n',path);
end
return
