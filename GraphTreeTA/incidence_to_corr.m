function [corr] = incidence_to_corr(M)
%function [corr] = incidence_to_corr(M)
%incidence_to_corr takes an edge-vertex incidence matrix M built on the
%relationship between epochs (single data) and pairs of epochs (two data
%points) and builds the corresponding coefficient matrix for the pairs.
%The output matrix, corr, is used to determine the covariance matrix for
%the set of pairs data. 
%Elena Baluyut
%2014-09-06

%Build Laplacian from incidence matrix**
L = M*M'; %edge Laplacian in terms of pair relationships

%ndat = length(L);

%Build Correlation Matrix
De = incidence_to_degree(M'); %Degree in terms of pair relationships (transpose of current DD) M'

corr = inv(sqrt(De))*L*inv(sqrt(De)); %theory reference: Merris, R. (1994) "Laplacian Matrices of Graphs: A Survey"


%
return

