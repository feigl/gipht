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

try
    status1 = 0;
    mphstart
catch ME
    fprintf(1,'COMSOL mph server is not started. Trying to restart...\n');
    switch computer
        case 'MACI64'
            [status1, output] = system('/Applications/COMSOL54/Multiphysics/bin/comsol mphserver &');
        otherwise
            error('Unknown computer');
    end
end

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
    warning('Could not start mphserver');
end


return
end




