function [year4s, months, days] = dyear2yyyymmdd(dyears)
%function [year4s, months, days] = dyear2yyyymmdd(dyears)
% given decimal year, return 4-digit year, month, day
% Should replace with function 'decyear' from Aerospace toolbox in Matlab release 2014a
% Examples
% [year4,month,day] = dyear2yyyymmdd(dyear(2016,3,3))
%   should return [2016, 3, 3]
% [year4,month,day] = dyear2yyyymmdd(dyear(2014,3,3))
%   should return [2014, 3, 3]
% datestr(datenum(year4,month,day))
%   03-Mar-2014
% 2014-06-19 Kurt Feigl

n = numel(dyears);
year4s = nan(n,1);
months = nan(n,1);
days   = nan(n,1);
for i=1:numel(dyears)
    dyear = dyears(i);
    year4 = floor(dyear);
    month = 1;
    day = 1;
    nday0 = datenum(dyear,1,1);
    nday1 = datenum(dyear,12,31);
    while nday1 - nday0 > 0
        [year4, month, day] = datevec(nday0);
        nday1 = datenum(year4,month,day);
        nday0 = nday0+1;
    end
    year4s(i) = year4;
    months(i) = month;
    days(i)   = day;
end

return
