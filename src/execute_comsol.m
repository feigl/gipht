function [model,solution_epochs] = execute_comsol(mphname)
%function solution_epochs = execute_comsol(mphname)
% run comsol mph file returning epochs of solution
% 20130726 Kurt Feigl

verbose = 1;
if verbose == 1;
    fprintf(1,'Loading %s\n',char(mphname));
end

model = mphload(char(mphname));

info0 = mphsolutioninfo(model);

if verbose == 1
    %info0
    fprintf(1,'Starting to calculate %s\n',char(mphname))
end;

tic1=tic;
model.sol('sol1').runAll;

% 20150901 call function
mphsave(model,char(mphname));

if verbose == 1
    fprintf(1,'Finished calculating %s in %.f seconds\n',char(mphname),toc(tic1));
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

