%% plot_comsol_optimization
% 20191111 Kurt Feigl 

% read the file
filename = 'ObjectiveTable6.csv'
% T=readtable(filename,'HeaderLines',6);
OPTS = detectImportOptions(filename);
OPTS.VariableNamesLine=5;
OPTS.VariableDescriptionsLine = 5;
OPTS.CommentStyle = {'%'};
OPTS.VariableNames= {'P1', 'Objective'}
T=readtable(filename,OPTS);
summary(T)
S=table2struct(T,'ToScalar',true)
figure
plot(S.P1/1.e6,S.Objective,'r*');
xlabel('P1 [MPa]');
ylabel('Objective function [dimless]');
title(strcat(filename);
printpdf(sprintf('%s',strrep(mfilename,'.m','.pdf')));
