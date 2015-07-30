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

% Find calendar dates
[yr mn dy] = dyear2date(tu);

% Find corresponding month string
mn_s = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% Print to string
for i = 1:numel(tu)
    cal_string{i} = sprintf('%4d-%s-%2d', yr(i), mn_s{mn(i)}, dy(i));
end

end

