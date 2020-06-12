function k=thermalcond(T)
% k=thermalcond(T)
%
% Warning: this function is just an example.
% The calculations doesn't reflect any real physical property.
%
% The input T is a (column) vector. 
% The output k must have exactly the same size as T.

k = 3-0.01*(T-298);
