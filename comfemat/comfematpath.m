function comfematpath
% set path for ComFEmat
home = getenv('HOME');
comfematv = 'comfemat';
comfemat_home = strcat(home,filesep,comfematv)
if exist(comfemat_home,'dir') == 7
    fprintf(1,'Setting environment variable GIPHT_HOME  to %s\n',comfemat_home);
    setenv('GIPHT_HOME',comfemat_home);
    fprintf(1,'Environment variable GIPHT_HOME is now set to %s\n',getenv('GIPHT_HOME'));
    p = strcat(...
         comfemat_home,filesep,'src',        pathsep ...
        );
    addpath(p,'-BEGIN');
    fprintf(1,'Matlab command search path is now: %s\n',path);
end
return
