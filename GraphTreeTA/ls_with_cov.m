function [x, sig, mse, Vx, Ginv] = ls_with_cov(A, B, V, Ginv)
% function [x, sig, mse, Vx] = ls_with_cov(A, B, V, Ginv)
% Gives the least-squares solution to A*x = B with the inverse data covariance matrix V as the weighting matrix.
%
% INPUTS:
%   A - design matrix, size[n - by - m]
%   B - data vector, size[n - by- 1]
%   V - data covariance matrix, size[n - by - n]
%   Ginv - optional inverse
% OUTPUTS:
%   x - least-squares estimated solution for unknown parameter vector, size
%       [m - by 1]
%   sig - vector of uncertainties for estimated solution
%   mse - mean squared error (or variance of unit weight) of estimated
%         solution
%   Vx - covariance matrix of estimated parameters
%
% Elena C. Baluyut, UW-Madison
% 2014-9-19
% 20200424 add optional arguments for Ginv

narginchk(3,4)
nargoutchk(3,5);

% calculate inverse matrices
if exist('Ginv','var') == 0 
    doinv = 1;
elseif numel(Ginv) < 1
    doinv = 1;
else
    [nGinv,mGinv] = size(Ginv);
    if nGinv == numel(B) 
        doinv = 0;
    else
        doinv =1;
    end
end

% calculate inverses
if doinv == 1
    % Find weighting matrix
    [S, rs] = pinveb(V); % calculates pseudoinverse similarly to pinv, but returns L-curve and rank
    % S = pinv(V);       % use when L-curve, etc. not needed
    
    % Calculate least-squares estimator
    [Ginv, rga] = pinveb(A'*S*A); % calculates pseudoinverse similarly to pinv, but returns L-curve and rank
    if rga > 10e6
        fprintf(1,'In %s Matrix is nearly singular (rga = %.4E), consider reducing the number of parameters\n,',mfilename,rga);
    end
    %  Ga = pinv(A'*S*A);        % use when L-curve, etc. not needed   
end

%size(Ginv)

x = Ginv*A'*S*B;

% Calculate statistics
[n, m] = size(A);
dfe = n-m; % degrees of freedom
r = B - A*x; % residuals
mse = (r'*S*r)./ dfe; %mean squared error estimate from formula 9.64 in strang and borre (1997) pg 344

% parameter covariance matrix
%Vx = pinv(A'*S*A)*mse; % no need to calculate inverse again
Vx = Ginv*mse; % 
% parameter uncertainty
%sig = sqrt(diag(pinv(A'*S*A)*mse)); % no need to calculate inverse again
sig = sqrt(diag(Vx)); 

return
end

