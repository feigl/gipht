function Y = heaviside1(X)
%HEAVISIDE1    Step function.
%    HEAVISIDE1(X) is 0 for X < 0, 1 for X >= 0
%    HEAVISIDE1(X) is not a function in the strict sense.
% WARNING: this function is differs from the Matlab built-in function
% HEAVISIDE in its behavior for X = 0
% Kurt Feigl 2011-OCT-20
Y = zeros(size(X));
Y(X >= 0) = 1;
return