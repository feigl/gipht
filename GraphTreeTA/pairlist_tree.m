function [ PN ] = pairlist_tree(pairlst, trees, rs)
%pairlist_tree: Using pairs, trees, and uncertainties, returns a matrix organizing the information for a minimum spanning tree.  Cycles and repeated pairs are removed according to Prim's algorithm.'
%   output matrix PN = [tree# pair# tm# ts# rsofpair#]

[ntrees, mdummy] = size(trees);
    PN = [];
for n = 1:ntrees
    iok = isfinite(trees(n,:));
    nodes = trees(n,iok); 
    for i = 1:numel(nodes)
        [i_tm i_tum] = find(pairlst(:,2) == nodes(i));
        %[i_ts i_tus] = find(pairlist(:,3) == nodes(i));
        if isempty(i_tm) == 1
            [i_tm i_tum] = find(pairlst(:,3) == nodes(i));
        end
        PN(end+1:end+numel(i_tm),1:5) = [n*ones(numel(i_tm),1) i_tm pairlst(i_tm,2) pairlst(i_tm, 3) rs(i_tm)];
%         for i = 1:numel(i_tm)
%         [i_cycles] = find(pairlist(i_tm, 2) == pairlist(i, 3));
        
        
    end

end

%Remove any repeated pairs 
[pp pq] = size(PN); RM = []; 
for l = 1:pp 
    [rmi] = find(PN(l, 2) == PN(:, 2));
    RM(end+1:end+numel(rmi)-1, :) = rmi(2:end);
end
PN(RM, :) = [];

return

