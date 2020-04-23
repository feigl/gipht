function [ pairs ] = pairlist( Q )
% function [ pairs ] = pairlist( Q )
% Given edge-vertex incidence matrix Q of a directed graph, finds the
% corresponding edges.
% For InSAR, the edges in the graph correspond to pairs (interferograms), 
% the vertices to epochs ("acquisition dates" for individual images)
% 
% INPUT:
%   Q     - edge-vertex inidence matrix 
% OUTPUT:
%   pairs - matrix containing information of pairs:
%           [index of pair, index of master, index of slave]
%   
% 2015-02-10 Elena C. Baluyut, UW-Madison
% 2020-03-26 Kurt Feigl clarify documentation.

% Initialize
[nedges nvertices] = size(Q); % get number of pairs
pairs = zeros(nedges, 3);  % set up storage matrix 

% Loop through pairs
for i = 1:nedges
    I = find(Q(i, :) ~= 0);     % find nonzero column indices 
    pairs(i,1) = i;             % store pair index - implied, not necessary
    pairs(i,2) = I(1);          % store index initial vertex of edge
    pairs(i,3) = I(2);          % store index terminal vertex of edge
end

return
end