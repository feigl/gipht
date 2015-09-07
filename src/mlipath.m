function mlipath
% add Matlab Livelink Interface to Matlab search path
mli_path = '/usr/local/comsol51/mli'
if exist(mli_path,'dir') == 7
    addpath(mli_path,'-BEGIN');
    fprintf(1,'Matlab command search path is now: %s\n',path);
end
return
