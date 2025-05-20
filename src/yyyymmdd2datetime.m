function [date_time] = yyyymmdd2datetime(yyyymmdd)
%function [date_time] = yyyymmdd2datetime(yyyymmdd)
% convert year (YYYYMMDD) to Matlab datetime



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

date_time = datetime(yyyy,mm,dd,'Format','yyyyMMdd','TimeZone','UTC');
% date_time.TimeZone='UTC';
% date_time.Format='yyyy-MM-dd''T''HH:mmXXX';


return

