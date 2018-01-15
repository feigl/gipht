function matlabdatetime = decyear2datetime(decyear)
%function dy = dyear(year, month, day)
% given decimal year, return Matlab datetime class
% Kurt Feigl 20180115
matlabdatetime = datetime(1, 1, 1) + years(decyear - 1);
return
