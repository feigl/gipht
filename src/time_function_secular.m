function ft = time_function_secular(tepochs, tref)
% return value of time function f(t)
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          trate  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch
%                    Heaviside step function
%          ft(t)    = t - tref
%
% 2011-OCT-17 Kurt Feigl
ft = zeros(size(tepochs));
itime=find(isfinite(tepochs));
ft(itime) = tepochs(itime)-tref;
return
end

