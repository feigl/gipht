function [corr] = incidence_to_corr(Q)
%function [corr] = incidence_to_corr(Q)
%incidence_to_corr takes an edge-vertex incidence matrix Q built on the
%relationship between epochs (single data) and pairs of epochs (two data
%points) and builds the corresponding coefficient matrix for the pairs.
%The output matrix, corr, is used to determine the covariance matrix for
%the set of pairs data.
%
% INPUT:
%   Q - edge-vertex incidence matrix (n (# pairs) - by - m (# epochs))
% OUTPUT:
%   corr - correlation matrix
%
%Elena C. Baluyut, UW-Madison
%2014-09-06

%Build Laplacian from incidence matrix**
L = Q*Q'; %edge Laplacian in terms of pair relationships 

%Build Correlation Matrix
De = incidence_to_degree(Q'); %Degree in terms of pair relationships (transpose of current DD) M'

corr = inv(sqrt(De))*L*inv(sqrt(De)); %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"

%
return

