function [ Gb, B ] = ber_mat( Q, tu)
% function [ output_args ] = ber_mat( Q, tu)
% Define Berardino et al. (2002) design matrix
%
% INPUTS:
%   Q -  edge-vertex incidence matrix
%   tu - vector of epochs
%
% OUTPUT:
%   Gb - Berardino design matrix
%
% Elena C. Baluyut, UW-Madison
% 2015-07-28

% Initialize
[ndat, mepochs] = size(Q);

% Build pair-vertex matrix
Del = zeros(mepochs-1, mepochs);
for j = 1:mepochs-1
    Del(j,j) = -1;  
    Del(j, j+1) = 1; 
end

% Build pair-rate incidence matrix 
B = Q*pinv(Del);

% Build diagonal time-interval matrix 
T = diag(tu(2:end)-tu(1:end-1));

% Define Berardino et al. (2002) design matrix
Gb = B*T; 
end

