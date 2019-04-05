function [caldate] = dyear2caldate(dyears)
%function [caldate] = dyear2caldate(dyears)
% given decimal year, return calendar date in YYYYMMDD
% 20180701 Elena C Reinisch

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
   caldate = year4s*1e4+months*1e2+days;
return
