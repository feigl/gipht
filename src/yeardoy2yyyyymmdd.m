function [year,month,day] = yeardoy2yyyyymmdd(yy,doy)
%function [year,month,day] = yeardoy2yyyyymmdd(yy,doy)
% convert year (YY) and day of year (DOY) to year, month, day

narginchk(2, 2);

if yy-floor(yy) > eps
    warning(sprintf('Truncating non-integer year (%f) to integer',yy));
    yy = floor(yy);
end

if yy < 0
    error(sprintf('negative year (%d)',yy));
end

if yy < 50
    yy = 2000 + yy;
elseif (50 <= yy) && (yy < 100)
    yy = 1900 + yy;
end

year = yy;

day1 = datetime(year,1,1) + days(doy-1);
month = day1.Month;
day = day1.Day;

return

