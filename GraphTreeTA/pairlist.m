function [ pairs ] = pairlist( Q )
% function [ pairs ] = pairlist( Q )
% Using the edge-vertex incidence matrix Q, finds the corresponding epochs for each pair.
%
% INPUT:
%   Q     - edge-vertex inidence matrix 
% OUTPUT:
%   pairs - matrix containing information of pairs:
%           [index of pair, index of master, index of slave]
%   
% Elena C. Baluyut, UW-Madison
% 2015-02-10

% Initialize
[ndat mdummy] = size(Q); % get number of pairs
pairs = zeros(ndat, 3);  % set up storage matrix 

% Loop through pairs
for i = 1:ndat
    I = find(Q(i, :) ~= 0);     % find nonzero column indices 
    pairs(i,1) = i;             % store pair index
    pairs(i, 2) = I(1);         % store index of master epoch of pair
    pairs(i,3) = I(2);          % store index of slave epoch of pair
end

return

