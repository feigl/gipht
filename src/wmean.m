function [d_mean,d_std,sqrt_mse,s_scaled] = wmean(d,s)
%function [d_mean,d_std,sqrt_mse,s_scaled] = wmean(d,s)
% given a column vector of data, return weighted mean and scaled sigmas
% calls Matlab routine LSCOV - see help lscov
% Example:
%
% sample a random variable with zero mean and unit standard deviation
% z=randn([10000,1]);
% Find its mean
% mean(z)
% ans =
%     0.0051

% [d_mean,d_std,sqrt_mse] = wmean(z,ones(size(z)))
% d_mean =
%     0.0051
% d_std =
%     0.0099
% sqrt_mse =
%     0.9854
% 
% [d_mean,d_std,sqrt_mse] = wmean(z,2.0*ones(size(z)))
% d_mean =
%     0.0051
% d_std =
%     0.0099
% sqrt_mse =
%     0.4927
% Now try including a missing value
% >> z(end+1)=NaN;
%
% Results should be the same
% >> [d_mean,d_std,sqrt_mse] = wmean(z,2.0*ones(size(z)))
% Warning: Omitting 1 values that are not finite
%  
% > In wmean at 124
% d_mean =
%     0.0051
% d_std =
%     0.0099
% sqrt_mse =
%     0.4927
%
%
% 2011-OCT-06 Kurt Feigl 

% 
% help lscov
%  LSCOV Least squares with known covariance.
%     X = LSCOV(A,B) returns the ordinary least squares solution to the
%     linear system of equations A*X = B, i.e., X is the N-by-1 vector that
%     minimizes the sum of squared errors (B - A*X)'*(B - A*X), where A is
%     M-by-N, and B is M-by-1.  B can also be an M-by-K matrix, and LSCOV
%     returns one solution for each column of B.  When rank(A) < N, LSCOV
%     sets the maximum possible number of elements of X to zero to obtain a
%     "basic solution".
%  
%     X = LSCOV(A,B,W), where W is a vector length M of real positive weights,
%     returns the weighted least squares solution to the linear system A*X =
%     B, i.e., X minimizes (B - A*X)'*diag(W)*(B - A*X).  W typically
%     contains either counts or inverse variances.
%  
%     X = LSCOV(A,B,V), where V is an M-by-M symmetric (or hermitian) positive
%     definite matrix, returns the generalized least squares solution to the
%     linear system A*X = B with covariance matrix proportional to V, i.e., X
%     minimizes (B - A*X)'*inv(V)*(B - A*X).
%  
%     More generally, if V is full, it can be positive semidefinite, and LSCOV
%     returns X that is a solution to the constrained minimization problem
%  
%        minimize E'*E subject to A*X + T*E = B
%          E,X
%  
%     where T*T' = V.  When V is semidefinite, this problem will have a
%     solution only if B is consistent with A and V (i.e., in the column
%     space of [A T]), otherwise LSCOV returns an error.
%  
%     By default, LSCOV computes the Cholesky decomposition of V and, in
%     effect, inverts that factor to transform the problem into ordinary
%     least squares.  However, if LSCOV determines that V is semidefinite, it
%     uses an orthogonal decomposition algorithm that avoids inverting V.
%  
%     X = LSCOV(A,B,V,ALG) allows you to explicitly choose the algorithm used
%     to compute X when V is a matrix.  LSCOV(A,B,V,'chol') uses the Cholesky
%     decomposition of V.  LSCOV(A,B,V,'orth') uses orthogonal decompositions,
%     and is more appropriate when V is ill-conditioned or singular, but is
%     computationally more expensive.  'orth' is not allowed when any of the
%     inputs are sparse.
%  
%     [X,STDX] = LSCOV(...) returns the estimated standard errors of X.  When
%     A is rank deficient, STDX contains zeros in the elements corresponding
%     to the necessarily zero elements of X.
%  
%     [X,STDX,MSE] = LSCOV(...) returns the mean squared error.
%  
%     [X,STDX,MSE,S] = LSCOV(...) returns the estimated covariance matrix
%     of X.  When A is rank deficient, S contains zeros in the rows and
%     columns corresponding to the necessarily zero elements of X.  LSCOV
%     cannot return S if it is called with multiple right-hand sides (i.e.,
%     size(B,2) > 1).
%  
%     The standard formulas for these quantities, when A and V are full rank,
%     are:
%  
%        X = inv(A'*inv(V)*A)*A'*inv(V)*B
%        MSE = B'*(inv(V) - inv(V)*A*inv(A'*inv(V)*A)*A'*inv(V))*B./(M-N)
%        S = inv(A'*inv(V)*A)*MSE
%        STDX = sqrt(diag(S))
%  
%     However, LSCOV uses methods that are faster and more stable, and are
%     applicable to rank deficient cases.
%  
%     LSCOV assumes that the covariance matrix of B is known only up to a
%     scale factor.  MSE is an estimate of that unknown scale factor, and
%     LSCOV scales the outputs S and STDX appropriately.  However, if V is
%     known to be exactly the covariance matrix of B, then that scaling is
%     unnecessary.  To get the appropriate estimates in this case, you should
%     rescale S and STDX by 1/MSE and sqrt(1/MSE), respectively.
%  
%     Class support for inputs A,B,V,W:
%        float: double, single
%  
%     See also mldivide, slash, lsqnonneg, qr.
% 
%     Reference page in Help browser
%        doc lscov

nargchk(2,2,nargin);
nargchk(1,4,nargout);

d_mean = NaN;
s_scaled = nan(size(d));
sqrt_mse=NaN;
d_std = NaN;

n = numel(d);
if numel(s) ~= n  
    error(sprintf('size mismatch: %d %d\n',numel(s),n));
end

% prune out NaN values
iok = find(isfinite(d));
iok = intersect(iok,find(isfinite(s)));
if numel(iok) < n
    warning(sprintf('Omitting %d values that are not finite\n',n-numel(iok)));
    d=d(iok);
    s=s(iok);
    n=numel(d);
end
if n > 1
% data values
B = reshape(d,n,1);
% weights are inverse variances
W = 1./(s.^2);
W = reshape(W,n,1);
% design matrix
A = ones(n,1);

% solve
[d_mean,d_std,mse] = lscov(A,B,W);

% scale factor
sqrt_mse = sqrt(mse);

% re-scale uncertainties
s_scaled = s * sqrt_mse;

end
return

end

