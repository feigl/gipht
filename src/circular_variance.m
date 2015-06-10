function S0 = circular_variance ( theta )
% circular_variance (theta)
% given angles theta in radians,
% return circular variance S0 in dimesionless units
% Kurt Feigl
% from Mardia 1972
%
xbar = mean_direction(theta);
n = numel(theta);
S0 = 1-sum(cos(theta-xbar))/n;
return
