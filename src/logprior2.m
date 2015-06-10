function lp=logprior2(params,lobounds,upbounds)
%function lp=logprior2(param,lobound,upbound)
% given parameter values, lower bounds and upper bounds, return log
% probability of model
% 2014
% from logprior.m
% Example 11.4
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber

% number of parameters
mparams = numel(params);

% increment count if in bounds
ngood = 0;
for i=1:mparams
    if params(i) >= lobounds(i) && params(i) <= upbounds(i)
        ngood = ngood+1;
    end
end

% if all in bounds, then the candidate is acceptable
if ngood == mparams
    lp=0;
else
    lp=-Inf;
end

return
end
