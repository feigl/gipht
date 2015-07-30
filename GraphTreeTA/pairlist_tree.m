function [ PN ] = pairlist_tree(pairlst, trees, rs)
% function [ PN ] = pairlist_tree(pairlst, trees, rs)
% Using pairs, components, and uncertainties, returns a matrix organizing the information for a minimum spanning tree.  
% Cycles and repeated pairs are removed according to Prim's algorithm.
%
% INPUTS:
%   trees - matrix of components in data set
%   rs    - uncertainty of pairwise data
%
% OUTPUT: 
%   PN    - matrix containing information of pairs: [tree# pair# tm# ts# rsofpair#]
%
% Elena C. Baluyut, UW-Madison
% 2015-02-09

% Initialize 
[ntrees, mdummy] = size(trees); % find number of trees in data set
PN = [];                        % set matrix for storing tree statistics 

% Loop through component
for n = 1:ntrees
    iok = isfinite(trees(n,:)); % find element locations of epoch indices in each component
    nodes = trees(n,iok);       % stores epoch index values in a new vector 
  
    % Loop through each epoch in component
    for i = 1:numel(nodes)
        [i_tm i_tum] = find(pairlst(:,2) == nodes(i)); % find other pairs where epoch is slave
        if isempty(i_tm) == 1
            [i_tm i_tum] = find(pairlst(:,3) == nodes(i)); % find other pairs where epoch is master
        end
        PN(end+1:end+numel(i_tm),1:5) = [n*ones(numel(i_tm),1), i_tm, pairlst(i_tm,2), pairlst(i_tm, 3), rs(i_tm)]; % store pairs sharing same epoch
    end
    
end

% Remove any repeated pairs 

[pp pq] = size(PN); RM = [];                     % initialize matrix for removed indexes
for l = 1:pp 
    [rmi] = find(PN(l, 2) == PN(:, 2));          % find pairs in cycles to be eliminated
    RM(end+1:end+numel(rmi)-1, :) = rmi(2:end);  % store pairs to be eliminated
end
PN(RM, :) = []; % remove unwanted pairs 

return

