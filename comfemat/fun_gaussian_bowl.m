function y = fun_gaussian_bowl(x,mu,sigma)
% return an upside-down Gaussian with mean mu and standard deviation, evaluate at x
% 2020/05/26 Kurt Feigl

y = normpdf(0,mu,sigma)-normpdf(x,mu,sigma);

return
end

