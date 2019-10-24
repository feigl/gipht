function [date_num] = yyyymmdd2datenum(yyyymmdd)
%function [date_num] = yyyymmdd2datenum(yyyymmdd)
% convert year (YYYYMMDD) to Matlab date number
% 20191024 Kurt Feigl
% examples:
% yyyymmdd2datenum('20191024')
% yyyymmdd2datenum(20191024)
% datenum(2019,10,24)



narginchk(1,1);

if ischar(yyyymmdd) ~= 1
    yyyymmdd = sprintf('%08d',yyyymmdd);
end

yyyy = str2double(yyyymmdd(1:4));
  mm = str2double(yyyymmdd(5:6));
  dd = str2double(yyyymmdd(7:8));
if yyyy-floor(yyyy) > eps
    warning(sprintf('Truncating non-integer year (%f) to integer',yyyy));
    yyyy = floor(yyyy);
end

if yyyy < 0
    error(sprintf('negative year (%d)',yyyy));
end

if yyyy < 50
    yyyy = 2000 + yyyy;
elseif (50 <= yyyy) && (yyyy < 100)
    yyyy = 1900 + yyyy;
end

date_num = datenum(yyyy,mm,dd,0,0,0);

return

