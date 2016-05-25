function dy = dyear(year, month, day)
%function dy = dyear(year, month, day)
%return decimal year
% Should replace with function 'decyear' from Aerospace toolbox in Matlab release 2014a
% Kurt Feigl
dy = years(datetime(year,month,day) - datetime(1, 1, 1)) + 1
return
