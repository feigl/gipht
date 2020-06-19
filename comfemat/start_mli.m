
function status1 = start_mli(verbose)
% start matlab live link with comsol on Mac
% 20200526 Kurt Feigl

narginchk(0,1);
nargoutchk(0,1);
if nargin < 1
    help(sprintf('%s',mfilename));
end

if nargin < 1
    verbose = 0;
end

%% Comsol server must be started
switch computer
    case 'MACI64'
        [status0, output] = system('pgrep mphserver')
        %if status0 == 0 && isempty(output) ==      
        pid = str2num(output);
        if status0 == 0 && numel(find(isfinite(pid) == true)) > 0
            fprintf(1,'COMSOL mph server is running.\n');
            status1 = 0;
        else
            fprintf(1,'COMSOL mph server is not started. Trying to start...\n');
            % 20200619 Must run in foreground. Do not add ampersand to this command.
            % 20200619 Name of option is different for Mac
            [status0, output] = system('/Applications/COMSOL54/Multiphysics/bin/comsol server')
            if status0 == 0 && isempty(output) == true
                fprintf(1,'COMSOL mph server successfully restarted.\n');
                status1 = 0;
            end
        end
    case 'GLNXA64'
        hostname = getenv('HOSTNAME');
        switch hostname
            case 'porotomo.geology.wisc.edu'
                fprintf(1,'Killing all comsol processes...\n');
                [status,output] = system('pkill -f comsol')
                fprintf(1,'Starting comsol server...\n');
                [status0, output] = system('module load comsol54; comsol server &; sleep 5')
                mlipath = '/usr/local/comsol54/mli';
                addpath(mlipath);
                mphstart;
                status1 = 0;
            otherwise
                hostname
                error('Unknown hostname');
        end
%         [status0, output] = system('pgrep -f mphserver');
%         %if status0 == 0 && isempty(output) == false
%         pid = str2num(output);
%         if status0 == 0 && numel(find(isfinite(pid) == true)) > 0
%             
%             fprintf(1,'COMSOL mph server is running.\n');
%             status1 = 0;
%         else
%             
%             [status,output] = system('pkill -signal 9 -f mphserver');
%             [status0, output] = system(sysstr);
%             addpath(mlipath);
%             mphstart;
%             status1 = 0;
%         end
%         

            otherwise
                hostname
                error('Unknown hostname');
        end
        [status0, output] = system('pgrep -f mphserver');
        %if status0 == 0 && isempty(output) == false
        pid = str2num(output);
        if status0 == 0 && numel(find(isfinite(pid) == true)) > 0
            
            fprintf(1,'COMSOL mph server is running.\n');
            status1 = 0;
        else
            fprintf(1,'COMSOL mph server is not started. Trying to start...\n');
            [status,output] = system('pkill -signal 9 -f mphserver');
            [status0, output] = system(sysstr);
            addpath(mlipath);
            mphstart;
            status1 = 0;
        end
        
    otherwise
        error('Unknown computer');
end

% % wait 3 seconds
% t = timer;
% t.StartDelay = 3;
% t.TimerFcn = @(myTimerObj, thisEvent)disp('3 seconds have elapsed');
% start(t);
% 
% 
% if status0 == 0
%     mphstart
% %     try
% %         status1 = 0;
% %         mphstart
% %     catch ME
% %         fprintf(1,'COMSOL mph server is not started. Trying to restart...\n');
% %         switch computer
% %             case 'MACI64'
% %                 [status1, output] = system('/Applications/COMSOL54/Multiphysics/bin/comsol mphserver &');
% %             case 'GLNXA64'
% %                 [status1, output] = system('module load comsol54;comsol mphserver &')
% %             otherwise
% %                 error('Unknown computer');
% %         end
% %     end
% end

if status1 == 0   
    %% import some utilities
    import com.comsol.model.util.*
    
    %% check on status
    status2 = which('mphload');
    if contains(status2,'mphload.p') == 1
        fprintf(1,'Matlab Livelink server is available.\n');       
        fprintf(1,'If this is new to you, consider the following commands:\n');
        fprintf(1,'  help mphload\n');
        fprintf(1,'  mphnavigator\n');
        fprintf(1,'  mphgetproperties\n');
        fprintf(1,'  mphdoc(\''mphload\'')\n');
        fprintf(1,'  model.param.varnames\n');
        fprintf(1,'  mphsolinfo(model,''soltag'',''sol1'')\n');
        fprintf(1,'  varNames = model.param.varnames\n');
        status1 = 0;
    else
        error('Matlab Livelink server is not available.\n');
    end
else 
    ME.message
    status1
    output
    fprintf(1,'Could not start mphserver. Consider following command:\n');
    fprintf(1,'system(''kill -9 `pgrep mphserver`'')');
    error('Could not start mphserver');
end


return
end




