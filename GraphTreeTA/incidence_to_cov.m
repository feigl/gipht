function [V] = incidence_to_cov(M, dsig, tfunc, ndat)
%function [V] = incidence_to_cov(M, dsig, tfunc, ntrees)
%incidence_to_cov takes an incidence matrix M built on the
%relationship between epochs (single data) and pairs of epochs (two data
%points) and uses the relationship between the Laplacian and incidence matrices
%of a graph to find the correlation matrix, corr.  Corr is then used to
%build the covariance matrix from error vector dsig
%M = incidence matrix, m x n, with  length(data) = m and length(pairs) = n
%dsig = data error vector
%tfunc = time function
%ntrees = number of constraints (only used for tfunc = 'pwl', else enter 0)
%Elena Baluyut
%2014-09-17

%Build Laplacian from incidence matrix**
L = M*M'; %edge Laplacian in terms of pair relationships

%ndat = length(L);

%Build Correlation Matrix
De = incidence_to_degree(M'); %Degree in terms of pair relationships (transpose of current DD) M'

corr = inv(sqrt(De))*L*inv(sqrt(De)); %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"

%Find Covariance Matrix
V = diag(dsig)*corr*diag(dsig);

%**alternatively, can also use L = De - A where De is degree matrix and A is adjacency matrix.  However,
%the current function for adjacency matrix finds only those corresponding to
%vertices (not the edge relationship needed for temporal adjustment) (- will
%be fixed in the future)
return


