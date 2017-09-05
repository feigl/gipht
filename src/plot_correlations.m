%% make scatter plots of parameters
% Kurt Feigl 20170830

%%find indices of free parameters
jfree = find(abs(p1-p0)>0);
varnames = strrep(pnames(jfree),'_',' ');
varnames = strrep(varnames,'Okada3','');
% varnames = extractAfter(pnames(jfree),'_');
% [nr,nc] = size(varnames)
% for i=1:nc
%     varnames_cell{i}=varnames(:,i);
% end


%% find indices of good solutions within confidence
iwithin = find(acosts1 < crit69);                %  threshold from test
%iwithin = find(acosts1 < quantile(acosts1,0.05)); % within 5 percent quantile

%% make the plots
close all


corrplot(trials(iwithin,jfree),'varNames',varnames);
%corrplot(trials(iwithin,jfree),'varNames',varnames_cell);
h=gcf;
set(h,'PaperType','a0');
%set(h,'PaperPositionMode','Manual');
% Position = get(h,'Position');
% set(h,'Position',10*Position);
print(h,'test.pdf','-dpdf','-r1200','-bestfit'); % print PDF 

%%
figure
gplotmatrix(trials(iwithin,jfree))