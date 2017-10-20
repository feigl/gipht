function ft = timefun_boxcar(tepochs, treference1, treference2)
% return value of boxcar time function f(t)
% f(t) = 0 for t < treference1
% f(t) = 1 for t >= reference1 and t < treference2
% f(t) = 0 for t >= treference2
% 20170915 Kurt Feigl

%    HEAVISIDE1(X) is 0 for X < 0, 1 for X >= 0
%    HEAVISIDE1(X) is neither continuous nor differentiable
% WARNING: this function is differs from the Matlab built-in function
% HEAVISIDE in its behavior for X = 0
if treference1 < treference2
   ft = heaviside1(tepochs - treference1) .* heaviside1(treference2-tepochs) * (treference2-treference1);
else
    treference1
    treference2
    error('Reference times out of order.\n');
end
return
end
