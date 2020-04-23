function [V] = incidence_to_cov2(Q, dsig)
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

[nrowsQ,ncolsQ] = size(Q);
if nrowsQ <= 1 || ncolsQ <= 1
    warning('Incidence matrix is too small.\n');
    V = nan;
    return  
end

[nrowsS,ncolsS] = size(dsig);
if nrowsS == 0 || ncolsS == 0
    warning('Dsig vector is too short.\n');
    V = nan;
    return  
end

if ncolsS > 1
    warning('Transposing Dsig vector.\n');
    dsig = dsig';
end



% Build Laplacian from incidence matrix**
L = Q*Q'; %edge Laplacian in terms of pair relationships - square
[nrowsL,ncolsL] = size(L);
if nrowsL ~= ncolsL
    nrowsL
    ncolsL
    error('Laplacian matrix is not square.\n');
end


% Build Correlation Matrix
De = incidence_to_degree(Q'); %Degree in terms of pair relationships (transpose of current DD) M'
[nrowsDe,ncolsDe] = size(De);
if nrowsDe ~= ncolsDe
    nrowsDe
    ncolsDe
    error('Correlation matrix is not square.\n');
end

%corr = inv(sqrt(De))*L*inv(sqrt(De)); %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"
DeInv = inv(sqrt(De));
% 2020/04/21 use pseudo inverse 
%DeInv = pinv(sqrt(De));
[nrowsDeInv,ncolsDeInv] = size(DeInv);
if nrowsDeInv ~= ncolsDeInv
    nrowsDeInv
    ncolsDeInv
    error('DeInv is not square.\n');
end

if nrowsL ~= ncolsDeInv || nrowsL ~= ncolsDeInv
    Q
    De
    L
    DeInv
    nrowsL
    nrowsDeInv
    ncolsL
    ncolsDeInv
    
    warning('matrices differ in size')
end

corr = DeInv * L * DeInv; %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"
%corr = DeInv * L * DeInv'; % try to guarantee symmetry?

% Find Covariance Matrix
V = diag(dsig)*corr*diag(dsig);
% 2020/04/21 try to guarantee symmetry
%V = diag(dsig)*corr*diag(dsig)';

return


