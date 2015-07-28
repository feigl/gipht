function [ pairs ] = pairlist( DD )
%pairlist: Using the incidence matrix DD, finds the corresponding epochs for each pair.
%   Detailed explanation goes here
[m ndummy] = size(DD);
pairs = zeros(m, 3);
for i = 1:m
    I = find(DD(i, :) ~= 0);
    pairs(i,1) = i;
    pairs(i, 2) = I(1);
    pairs(i,3) = I(2);
end

