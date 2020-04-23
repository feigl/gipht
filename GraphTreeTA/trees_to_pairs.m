function ipairs = trees_to_pairs(itrees,iedges)
% given a list of trees in a directed graph, generate a list indices to pairs
%     itrees == list of trees, giving list of indices to vertices in each tree, 
%               one row per tree
%     iedges == (n,2) array listing vertices, starting and ending
% outputs: 
%     ipairs == list of pairs
% Reference:
% Reinisch, E.C., Cardiff, M. & Feigl, K.L. Graph theory for analyzing
% pair-wise data: application to geophysical model parameters estimated
% from interferometric synthetic aperture radar data at Okmok volcano,
% Alaska. J Geod 91, 9?24 (2017).
% https://doi.org/10.1007/s00190-016-0934-5
%
% 2012003228 Kurt Feigl

[ntrees,nvertices] = size(itrees);
[nedges,ncols] = size(iedges);

if ncols ~= 2
    error('Number of columns in list of edges must equal 2');
end


%fprintf(1,'The following lines should reproduce the edge list.\n')
kount = 0;
ipairs=nan(nedges,1);
for i=1:nedges
    % list of vertices in this tree
    j1 = iedges(i,1);
    j2 = iedges(i,2);
    for j=1:ntrees
        itree1 = itrees(j,:);
        % strip out nan values
        itree1 = itree1(isfinite(itree1));
        k1 = find(itree1 == j1);
        k2 = find(itree1 == j2);
        if numel(k1) == 1 && numel(k2) == 1         
%             fprintf(1,'%d ',itree1(k1));
%             fprintf(1,'%d ',itree1(k2));
%             fprintf(1,'\n');          
            kount = kount+1;
            ipairs(kount) = i;
        end
        ipairs=colvec(ipairs(isfinite(ipairs)));
    end
end
return
end

