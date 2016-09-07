function ft = timefun_step(tepochs, treference)
% return value of step time function f(t)
ft = heaviside1(tepochs - treference);
end
