function [ i_rept, nc ] = MST(trees, n, sub_PN )
% function [ i_rept, nc ] = MST(trees, n, sub_PN )
% Finds indices that should be eliminated to form the minimum spanning tree (Prim's
% algorithm).
%
% INPUTS:
%   trees  - matrix of components in data set
%   n      - index of distinct component 
%   sub_PN - submatrix of PN with only data for the nth component (see
%            'MSF.m' and 'pairlist_tree.m')
%
% OUTPUTS:
%   i_rept - indices of repeated pairs that are not to be included in MST
%   nc     - number of cycles in the nth component
%
% Elena C. Baluyut, UW-Madison
% 2015-02-19

% Initialize
iok = isfinite(trees(n,:)); % find element locations of epoch indices in each component
nodes = trees(n,iok); % store epoch index values in a new vector
ikeep = []; % initialize storage for indices to keep
nc = 0; %number of cycles

% loop through epoch indices in component 
for i = 1:numel(nodes)
    pairs =  sub_PN(find(sub_PN(:, 3) == nodes(i)),:); %find pairs with common master
    if isempty(pairs) == 1
        pairs = sub_PN(find(sub_PN(:,4) == nodes(i)), :); %find pairs with common slave
    end
    for k = 1:numel(pairs(:,4))
        [ex ex_v] = find(sub_PN(:,3) == pairs(k, 4)); %find slaves of node(i) that are also masters
        if isempty(ex) == 1
            % ikeep(end+1, :) = pairs(k, 2);
        else
            slaves_ex = sub_PN(ex, 4);  % find slave epochs corresponding to indices that are both master and slave
            for y = 1:numel(pairs(:,4)) % loop through slaves with shared epoch
                for h = 1:numel(ex)
                    rep = find(pairs(y,4) == slaves_ex(h)); % find indices of pairs that form cycle
                    if isempty(rep) == 0
                        % build matrix that stores cycle
                        cycle =   [pairs(y,2) pairs(y, 3) pairs(y, 4) pairs(y, 5);...
                            pairs(k, 2) pairs(k, 3) pairs(k, 4) pairs(k, 5);...
                            sub_PN(ex(h), 2) sub_PN(ex(h), 3) sub_PN(ex(h), 4), sub_PN(ex(h), 5)];
                        [index order] = sort(cycle(:,end)); % sort cycle based on uncertainty
                        ikeep(end+1,:) = cycle(order(3),1); % store index of pair with highest uncertainty (to be removed)
                    end
                    nc = nc + numel(rep); % count and store number of cycles per component 
                end
            end
        end
    end
end
[b,m1,n1] = unique(ikeep,'first'); % store index of first occurence of index to be removed (prevents counting index more than once)
[c1,d1] =sort(m1); % find indices of ikeep(m1)
i_rept = b(d1); % store repeated indices

end
