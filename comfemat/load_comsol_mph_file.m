function [model,info0] = load_comsol_mph_file(fileNameMPH,verbose)
%function [model,info0,info1] = load_comsol_mph_file(fileNameMPH,verbose)
% set up MATLAB livelink with Comsol and load an MPH file
% If this is new to you, consider the following commands:
% help mphload
% mphnavigator
% mphgetproperties
% mphdoc('mphload')
% model.param.varnames
% mphsolinfo(model,'soltag','sol1')
% varNames = model.param.varnames    

% 20200420 Kurt Feigl

narginchk(0,2);
nargoutchk(0,2);
if nargin < 1
    help(sprintf('%s'\n',mfilename));
end

if nargin < 2
    verbose = 0;
end
% 
% 
% %% Comsol server must be started
% 
% try
%     status = 0;
%     mphstart
% catch ME
%     fprintf(1,'COMSOL mph server is not started. Trying to restart...\n');
%     switch computer
%         case 'MACI64'
%             [status, output] = system('/Applications/COMSOL54/Multiphysics/bin/comsol mphserver &');
%         otherwise
%             error('Unknown computer');
%     end
% end
% if status ~= 0
%     ME.message
%     status
%     output
%     error('Could not start mphserver');
% end
% 
% %% import some utilities
% import com.comsol.model.util.*
% 
% %% check on status
% status = which('mphload');
% if contains(status,'mphload.p') == 1
%     fprintf(1,'Matlab Livelink server is available.\n');
%     
%     fprintf(1,'If this is new to you, consider the following commands:\n');
%     fprintf(1,'  help mphload\n');
%     fprintf(1,'  mphnavigator\n');
%     fprintf(1,'  mphgetproperties\n');
%     fprintf(1,'  mphdoc(\''mphload\'')\n');
%     fprintf(1,'  model.param.varnames\n');
%     fprintf(1,'  mphsolinfo(model,''soltag'',''sol1'')\n');
%     fprintf(1,'  varNames = model.param.varnames\n');    
% else
%     error('Matlab Livelink server is not available.\n');
% end

if nargin > 0
    model=mphload(fileNameMPH);
    
    % get some information
    info0 = mphsolinfo(model)
    
    
%     if verbose == 1
%         % get the names of the variables]
%         fprintf(1,'\nId Parameter_Value     [units]\n')
%         varNames = model.param.varnames;
%         for i=1:numel(varNames)
%             fprintf(1,'%3d %-32s %10.4G [%s]\n',i,varNames(i) ...
%                ,model.param.evaluate(sprintf('%s',varNames(i))) ...
%                ,model.param.evaluateUnit((sprintf('%s',varNames(i)))));
%         end
%         
%         %% list the names of the fields
%         fprintf(1,'\nId fieldName\n');
%         fieldTags = model.field.tags;
%         for i=1:numel(fieldTags)
%             data1 = mphgetfield(model.field(sprintf('%s',fieldTags(i))));
%             fieldNames{i} = sprintf('%s',char(data1.field));
%             fprintf(1,'%3d %-32s\n',i,fieldNames{i});
%         end
%     end
    
end


return
end



