function [ L ] = tikh_l_matrix( mparams, order )
% function [ L ] = tikh_l_matrix( mparams, order )
% Given the number of parameters in the model and the order for regularization,
% returns the Tikhonov roughening matrix. Supports zero-, first-, and second- order Tikhonov 
%
% INPUTS:
%   mparams - number of unknown parameters in linear system of equations
%   order   - order of regularization (0, 1, or 2)
%
% OUTPUT:
%   L       - Tikhonov roughening matrix
%
% Elena C. Baluyut, UW-Madison 
% 2015-07-01

if order == 0
    L = eye(mparams);
elseif order == 1
    L = zeros(mparams-1, mparams); % set size of L matrix
    for i = 1:mparams-1
        L(i, i) = -1; % assign -1 elements
        L(i, i+1) = 1; % assign +1 elements
    end
elseif order == 2
    L = zeros(mparams-2, mparams); % set size of L matrix
    for i = 1:mparams-2
        L(i, i) = 1; % assign second +1 elements
        L(i, i+1) = -2; % assign -2 elements
        L(i, i+2) = 1; % assign second +1 elements
    end
else 
    error('Function only supports zero-, first-, and second- order Tikhonov.  Please choose an appropriate order')
end

end

