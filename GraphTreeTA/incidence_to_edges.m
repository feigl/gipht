function [edges] = incidence_to_edges(Qiev)
% function [edges] = incidence_to_edges(Qiev)
% given incidence matrix of directed graph, find list of edges
% input:
%     Qiev  == Edge-vertex incidence matrix
% output:
%     edges == list of indices to edges
% Reference
% Reinisch, E.C., Cardiff, M. & Feigl, K.L. Graph theory for analyzing
% pair-wise data: application to geophysical model parameters estimated
% from interferometric synthetic aperture radar data at Okmok volcano,
% Alaska. J Geod 91, 9?24 (2017). https://doi.org/10.1007/s00190-016-0934-5
%
% 20200326 Kurt Feigl

[nedges, nvertices] = size(Qiev);
id0 = zeros(nedges,1); 
id1 = zeros(nedges,1); 
for i=1:nedges      
   ddcol = Qiev(i,:); 
   j=find(abs(ddcol)>0);
   id0(i) = min(j);  % index to first vertex
   id1(i) = max(j);  % index to second epoch
end
edges = [id0, id1];
return
end

