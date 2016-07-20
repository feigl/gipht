function [x, sig, mse, Vx] = ls_with_cov_bayes(m0, Vm, A, B, V)
% function [x, sig, mse, Vx] = ls_with_cov(A, B, V)
% Gives the least-squares solution to A*x = B with the inverse data covariance matrix V as the weighting matrix.
%
% INPUTS:
%   m0 - prior information on model
%   Vm - prior model covariance
%   A - design matrix, size[n - by - m]
%   B - data vector, size[n - by- 1]
%   V - data covariance matrix, size[n - by - n]
%
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

% Find weighting matrix    
 [S] = pinv(V); % calculates pseudoinverse similarly to pinv, but returns L-curve and rank
 
 
% Calculate least-squares estimator
 %[Ga, rga] = pinveb(A'*S*A); % calculates pseudoinverse similarly to pinv, but returns L-curve and rank
%   if rga > 10e6
%      warning('nearly singular, consider reducing the number of parameters')
%   end

    x = m0 + Vm*A'*pinv(V+A*Vm*A')*(B-A*m0); % formula from Tarantola & Nercessian 1984
  
% Calculate statistics
    [n, m] = size(A);
    dfe = n-m; % degrees of freedom
    r = B - A*x; % residuals
    mse = (r'*S*r)./ dfe; %mean squared error estimate from formula 9.64 in strang and borre (1997) pg 344
    Vx = pinv(A'*S*A)*mse; % parameter covariance matrix
    sig = sqrt(diag(pinv(A'*S*A)*mse)); % parameter uncertainty
                         
return

