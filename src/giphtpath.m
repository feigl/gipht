function giphtpath
%% set path for GIPHT
% 20200611 include comfemat

% switch computer
%     case 'GLNXA64'
%         hostname = getenv('HOSTNAME');
%         switch hostname
%             case 'porotomo.geology.wisc.edu'
%                 home = '/usr1/feigl';
%                 %home = '/usr1/ebaluyut';              
%             otherwise
%                 home = getenv('HOME');
%                 error('Unknown hostname');
%         end
%     case 'MACI64'       
%         home = '/Users/feigl';
%     otherwise
%         error('Unknown computer');
% end


% giphtv = 'gipht';
% gipht_home = strcat(home,filesep,giphtv)
% %gipht_home = '/Users/Elena/InSAR_development/TEST_GIPHT/gipht_porotomo20161026';
% if exist(gipht_home,'dir') == 7
%     fprintf(1,'Setting environment variable GIPHT_HOME  to %s\n',gipht_home);
%     setenv('GIPHT_HOME',gipht_home);
%     fprintf(1,'Environment variable GIPHT_HOME is now set to %s\n',getenv('GIPHT_HOME'));
%     p = strcat(...
%          gipht_home,filesep,'src',        pathsep ...
%         ,gipht_home,filesep,'pha2qls',    pathsep ...
%         ,gipht_home,filesep,'GraphTreeTA',pathsep ...
%         ,gipht_home,filesep,'extern',     pathsep ...
%         ,gipht_home,filesep,'utils',      pathsep ...
%         ,gipht_home,filesep,'comfemat',      pathsep ...
%         );
%     addpath(p,'-BEGIN');
%     fprintf(1,'Matlab command search path is now: %s\n',path);
% end

%restoredefaultpath; % factory default
home = getenv('HOME');
%https://github.com/feigl/gipht
%home = '/Volumes/GoogleDrive/Shared drives/Coso/gipht-master'
if numel(home) > 0
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'GraphTreeTA')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'extern')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'utils')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'src')),'-begin');
    %addpath(genpath('/Volumes/GoogleDrive/Shared drives/PoroTomo/Software/MatlabPoroTomo/UnitConversion_rev2'),'-begin');
else
    addpath(genpath('~feigl/gipht/GraphTreeTA'),'-begin');
    addpath(genpath('~feigl/gipht/extern'),'-begin');
    addpath(genpath('~feigl/gipht/utils'),'-begin');
    addpath(genpath('~feigl/gipht/src'),'-begin');
end

return
