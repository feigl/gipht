function [ i_msf, i_rep, Nc ] = MSF(ntrees, trees, PN)
% function [ i_msf, i_rep, Nc ] = MSF(ntrees, trees, PN)
% Finds minimum spanning forest (MSF) of data set based on Prim's Algorithm
% 
% INPUTS:
%   ntrees - number of components in data set
%   trees  - matrix of components in data set
%   PN     - matrix containing information of pairs (see function
%            'pairlist_tree.m')
%
% OUTPUTS:
%   i_msf - indices of epochs in the minimum spanning forest
%   i_rep - repeated indices of epochs to be eliminated
%   Nc    - vector the length of ntrees showing the number of cycles per each
%           component
%
% Elena C. Baluyut, UW-Madison
% 2015-02-19

% Initialize storage vectors
i_rep = [];
Nc = zeros(ntrees, 1);
[i_max mdummy] = size(PN); % find maximum number of pairs 
i_tot = [1:1:i_max]';   % find pair indices

% loop through components
for k = 1:ntrees
    i_start = find(PN(:,1) == k, 1, 'first'); % find start of list for component k
    i_last = find(PN(:,1) == k, 1, 'last'); % find end of list for component k
    sub_PN = PN(i_start:i_last, :); % select section of PN partaining to component k
    [ i_repmst, nc ] = MST(trees, k, sub_PN ); % find indices not included in minimum spanning tree of component k
    i_rep(end+1:end+numel(i_repmst), :) = i_repmst; % store indices not included in MST of k 
    Nc(k, 1) = nc; % store number of cycles in component k
end

i_rep = sort(i_rep); 
i_tot(i_rep) = []; % remove repeated indices from index list
i_msf = i_tot; % store indices of MST of all components for minimum spannig forest (MSF)

return


