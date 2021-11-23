function [s_hat, xi, beta, C_post] = lin_geostat_inv(y,X,R,Q,H)

%lin_geostat_inv.m: Function used to compute the best estimate of
%parameters given data for a linear forward model.
%   [s_hat, xi, beta,[s_postcov]] = lin_invert_solve(y,X,R,Q,H)
%   where:
%   OUTPUTS:
%       -s_hat is the best estimate of the parameters (m x 1)
%       -xi is the set of coefficients for the "basis functions" as
%       described in Snodgrass and Kitanidis 
%       These are multipliers for the matrix Q*H')
%       -beta is the set of drift function coefficients (k x 1)
%       -C_post is an m x m matrix containing the estimated posterior covariance matrix
%   INPUTS:
%       -y is the data (n x 1). If solving a nonlinear case, y should be
%       passed as y-h(s_tilde) + H*s_tilde
%       -X is the assignment or drift matrix (m x k)
%       -R is the data error covariance matrix (n x n)
%       -Q is the spatial covariance matrix of the parameters (n x n)
%       -H is the forward model design matrix (n x m) 
%       "H is an n by m matrix" Kitandis [2007] page 22

%
% Code by Michael Cardiff, 2009-2016
% 2021/09/11 Kurt Feigl - modify to check matrices first

% References:
% Kitanidis, P. K. (2007), 
% On Stochastic Inverse Modeling, in Subsurface Hydrology: Data
% Integration for Properties and Processes Geophysical Monograph Series 171, edited.
% https://agupubs-onlinelibrary-wiley-com.ezproxy.library.wisc.edu/doi/pdf/10.1029/171GM04 

% Cardiff, M. (2020), Bayesian Approaches to Inverse Problems, edited. 

% 1/2*(y-H*s)'*inv(R)*(y-H*s) + 1/2*(s-X*beta)'*inv(Q)*(s-X*beta)
% s_hat = X*beta + Q*H'*xi

num_reqout = 3;

%m = size(y,1);
ndata = size(y,1);  % number of data
mparam = size(X,1); % number of parameters
kdrift = size(X,2); % number of unknown "drift" parameters

% check dimensions
[nH,mH] = size(H);
if nH ~= ndata || mH ~= mparam
    error(sprintf('Design matrix H has incorrect dimensions %d %d\n',nH,mH));
end
% check dimensions of model covariance matrix
[nQ,mQ] = size(Q);
if nQ ~= mparam || mQ ~= mparam
    error(sprintf('Model covariance matrix has incorrect dimensions %d %d\n',nQ,mQ));
end
condQ=cond(Q)

% check dimensions of data covariance matrix
[nR,mR] = size(R);
if nR ~= ndata|| mR ~= ndata
    error(sprintf('Data covariance matrix has incorrect dimensions %d %d\n',nR,mR));
end
condR=cond(R)

% upper left hand corner of expanded design matrix in equation (9)
PSI = H*Q*H' + R;        

% off-diagonal term of of expanded design matrix in equation (9)
PHI = H*X;

%  \   Backslash or left division.
%       A\B is the matrix division of A into B, which is roughly the
%       same as INV(A)*B , except it is computed in a different way.
%       If A is an N-by-N matrix and B is a column vector with N
%       components, or a matrix with several such columns, then
%       X = A\B is the solution to the equation A*X = B. A warning
%       message is printed if A is badly scaled or nearly 
%       singular.  A\EYE(SIZE(A)) produces the inverse of A.

% This is equation (9) of Kitanidis, P. K. (2007), 
% On Stochastic Inverse Modeling, in Subsurface Hydrology: Data
% Integration for Properties and Processes Geophysical Monograph Series 171, edited.
% https://agupubs-onlinelibrary-wiley-com.ezproxy.library.wisc.edu/doi/pdf/10.1029/171GM04 

% xi_beta_vector = ([PSI, PHI; PHI', zeros(p,p)]) \...
% [y; zeros(p,1)];


% Left hand side of equation (9) of Kitanidis, P. K. (2007)
A9 = ([PSI, PHI; PHI', zeros(kdrift,kdrift)]);
condA9 = cond(A9)

% Right hand side of equation (9) of Kitanidis, P. K. (2007)
b = [y; zeros(kdrift,1)];

% least squares solution to equation (9) of Kitanidis, P. K. (2007)
%xi_beta_vector = A9 \ b;
xi_beta_vector = pinv(A9)* b;

% partition solution vector
xi = xi_beta_vector(1:ndata);

% estimate of drift parameters
beta = xi_beta_vector(ndata+1:ndata+kdrift);

% estimate of model parameters by Equation (9) of Kitanidis [2007]
s_hat = X*beta + Q*H'*xi;

% posterior model covariance matrix - original code
%s_postcov = Q - [H*Q; X']'*(A\[H*Q; X']);

% Equation (44) of Cardiff [2020]
[nHX,mHX] = size(H*X)
A44 = [(H*Q*H'+R), (H*X); (H*X)', zeros(kdrift)];
condA44=cond(A44)
% 2nd term of Equation (44) of Cardiff [2020]
s44 = ([H*Q; X'])' *pinv(A44) * [H*Q; X'];
[ns44,ms44] = size(s44);
% sum
C_post = Q - s44;

%To generate realizations:
% s_relz = s_hat + chol(s_postcov)'*randn(n,1);
%Similarly to:
% s_relz = mean + sqrt(var)*randn(n,1);

end