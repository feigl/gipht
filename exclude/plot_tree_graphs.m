function ierr = plot_tree_graphs(Tpair,Qiev,Edges,Weights,itrees,iedges)
nf=0; % number of figures
ierr = 0;
%% make title for plots
%titlestr = sprintf('%s %s %s',strtrim(timeSpan1),strtrim(dataSet),strtrim(criterion1));
titlestr = sprintf('%s',mfilename);

tu_cal   = unique([Tpair.mast;  Tpair.slav]);  % unique epochs
%% estimate the Yvalues for the graph from the weights
% this centers the yvalues
YvalEst = invert_incidence_matrix(Weights,Qiev,itrees);

%% make a graph using Matlab functions
DGRAPH = digraph(iedges(:,1),iedges(:,2),Weights,tu_str);
Nodes = table2array(DGRAPH.Nodes);

%% calculate importance of nodes
importance = centrality(DGRAPH,'authorities','Importance',abs(DGRAPH.Edges.Weight));
fprintf(1,'i,Nodes{i},tu_cal(i),importance(i) \n');
for i=1:nepochs
    fprintf(1,'%3d %8s %8d %10.4f\n',i,Nodes{i},tu_cal(i),importance(i));
end

%% is it a DAG?
idagness = isdag(DGRAPH);

%% plot trees
xlab = 'Date';
xunits = 'year';
xlabstr = sprintf('%s [%s]',xlab,xunits);
ylabstr = sprintf('%s [%s]',ylab,yunits);


% plot graph using Matlab function with labels
%'auto', 'circle', 'force', 'layered', 'subspace', 'subspace3', 'force3'
nf=nf+1;figure;
%layout = 'force';
layout = 'auto';
%layout = 'circle';
%layout = 'layered';
%layout = 'subspace';
plot(DGRAPH,'Layout',layout,'EdgeLabel',DGRAPH.Edges.Weight);
title(titlestr);
%                 printpdf(sprintf('%s_trees1_%s_%s_%s.pdf',mfilename ...
%                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
fname2 = 'digraph1';
PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
printpdf(PDFfileName);

% plot graph using Matlab function with labels
nf=nf+1;figure;
plot(DGRAPH,'EdgeLabel',DGRAPH.Edges.Weight,'XData',tu_dyear, 'YData',YvalEst);
xlabel(xlabstr);
ylabel(ylabstr);
title(titlestr);
%                 printpdf(sprintf('%s_trees2_%s_%s_%s.pdf',mfilename ...
%                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
fname2 = 'digraph2';
PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
printpdf(PDFfileName);



% plot graph using function from GraphTreeTA
nf=nf+1;Hfig = plot_trees3(tu_dyear, YvalEst, Qiev, itrees, xlab, ylabstr, titlestr);
%                 printpdf(sprintf('%s_trees3_%s_%s_%s.pdf',mfilename ...
%                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
fname2 = 'plot_trees3';
PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
printpdf(PDFfileName);

return
end

