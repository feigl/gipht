function [model,solution_epochs] = execute_comsol(mphname)
%function solution_epochs = execute_comsol(mphname)
% run comsol mph file returning epochs of solution
% 20130726 Kurt Feigl

verbose = 1;
if verbose == 1;
    fprintf(1,'Loading %s\n',mphname);
end

model = mphload(mphname);

info0 = mphsolutioninfo(model);

if verbose == 1
    %info0
    fprintf(1,'Starting to calculate %s\n',mphname)
end;

tic1=tic;
model.sol('sol1').runAll;
%     %model.sol('sol2').runAll; % not needed if Darcy and Poro solutions are combined
%    model.sol('sol1').updateSolution; % no need to remesh

%mphsave(model,sprintf('%s.mph',mphname));
%mphsave(model,mphname);

% 20140829 need to instance it 
model.save(mphname);

if verbose == 1
    fprintf(1,'Finished calculating %s in %.f seconds\n',mphname,toc(tic1));
end

info1 = mphsolinfo(model);
info2 = mphsolutioninfo(model);

% if verbose == 1
%     info1
%     info2
% end

% Comsol evaluation times in seconds
solution_epochs = info1.solvals;
return
end

