function dy = dyear(year, month, day)
%function dy = dyear(year, month, day)
%return decimal year
% Should replace with function 'decyear' from Aerospace toolbox in Matlab release 2014a
% Kurt Feigl
[doy,frac]=date2doy(datenum([year month day]));
dy = year+frac;
return
