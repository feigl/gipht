function [De] = incidence_to_degree(Q)
% function [De] = incidence_to_degree(Q)
% Given edge-vertex incidence matrix Q, returns degree matrix De (calculated column-wise).
% For ndat = (number of pairs) and me = (number of epochs)...
%   If Q is edge-vertex incidence matrix (ndat x me): 
%       vertex degree matrix = incidence_to_degree(Q)
%       edge degree matrix = incidence_to_degree(Q')
%   In GraphTeeTA, the edge degree matrix is used
%
% INPUT:
%   Q - edge-vertex incidence matrix (n (# pairs) - by - m (# epochs))
% OUTPUT:
%   De - degree matrix
%
% Elena C. Baluyut, UW-Madison
% 2014-08-27

De = diag(sum(abs(Q)));

return

