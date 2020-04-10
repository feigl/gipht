function decyear = datetime2dyear(caldate)
%function decyear = datetime2dyear(caldate)
%return decimal year from datetime object
% Elena Reinisch 20180329
decyear = years(caldate - datetime(1, 1, 1)) + 1;
return
