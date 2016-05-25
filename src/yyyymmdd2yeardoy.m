function [year,doy] = yyyymmdd2yeardoy(yy,month,day)
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

T0 = datetime(yy,1,1);
T1 = datetime(yy,month,day);
doy = days(T1 - T0) + 1;

return

