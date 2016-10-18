function [ cal_string ] = print_caldate( tu )
%function [ cal_string ] = print_caldate( tu )
%   Given a vector of epochs, returns a string vector with corresponding
%   calendar dates output as "YYYY-MON-DD"
% INPUT:
%   tu - vector of epochs
% OUTPUT:
%   cal_string - string vector containing calendar date as "YYYY-MON-DD"
%
% Elena C. Baluyut, UW-Madsion
% 2015 - 07 - 20

% check if dates are in calendar year format already
if round(tu(1)) ~= tu(1) 
    % if not find calendar dates
    [yr mn dy] = dyear2date(tu);
else % if in calendar dates
    tu_tmp = num2str(tu);
    yr = str2num(tu_tmp(:,1:4));
    mn = str2num(tu_tmp(:,5:6));
    dy = str2num(tu_tmp(:,7:8));
end

% Find corresponding month string
mn_s = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% Print to string
for i = 1:numel(tu)
    cal_string{i} = sprintf('%4d-%s-%2d', yr(i), mn_s{mn(i)}, dy(i));
end

end

