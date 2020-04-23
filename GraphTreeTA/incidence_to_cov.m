function [V] = incidence_to_cov(Q, dsig)
% function [V] = incidence_to_cov(Q, dsig)
% Takes an edge-vertex incidence matrix Q built on the
% relationship between epochs (single data) and pairs of epochs (two data
% points) and uses the relationship between the Laplacian and incidence matrices
% of a graph to find the correlation matrix, corr.  Corr is then used to
% build the covariance matrix from error vector dsig
%
% INPUTS:
%   Q    - edge-vertex incidence matrix, n x m, with  number of pairs = n and number of
%          epochs = m
%   dsig - data error vector
%
% OUTPUTS:
%   V    - covariance of pairwise data
%
% Elena C. Baluyut, UW-Madison
% 2014-09-17


% Build Laplacian from incidence matrix**
L = Q*Q'; %edge Laplacian in terms of pair relationships - square


% Build Correlation Matrix
De = incidence_to_degree(Q'); %Degree in terms of pair relationships (transpose of current DD) M'

corr = inv(sqrt(De))*L*inv(sqrt(De)); %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"


% Find Covariance Matrix
V = diag(dsig)*corr*diag(dsig);

return
end


