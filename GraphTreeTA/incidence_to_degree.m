function [De] = incidence_to_degree(DD)
%given incidence matrix DD, returns degree matrix De.
%For ndat = (number of pairs) and me = (number of epochs)...
%If DD is vertex-edge incidence matrix (me x ndat):
%   edge degree matrix = incidence_to_degree(DD)
%   vertex degree matrix = incidence_to_degree(DD')
% If DD is edge-vertex incidence matrix* (ndat x me): 
%   vertex degree matrix = incidence_to_degree(DD)
%   edge degree matrix = incidence_to_degree(DD')
%
%* as is the case in temporal adjustment
%
%Elena Baluyut 2014-08-27

De = diag(sum(abs(DD)));

return

