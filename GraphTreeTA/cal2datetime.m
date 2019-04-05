function [ datetime_vec ] = cal2datetime( tu )
% [ datetime_vec ] = cal2datetime( tu )
%  takes dates in calendar yyyymmdd and converts to datetime
%  Elena C Reinisch 20170503

 test = num2str(tu);
datetime_vec = datetime(str2num(test(:,1:4)), str2num(test(:, 5:6)), str2num(test(:, 7:8)));
return

