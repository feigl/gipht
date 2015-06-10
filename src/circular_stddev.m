function s0 = circular_stddev ( theta )
% circular_stddev (theta)
% given angles theta in radians,
% return circular standard deviation in radians
% Kurt Feigl
% from Mardia 1972
%
S0 = circular_variance(theta);
s0 = sqrt(-2*log(1-S0));
return
