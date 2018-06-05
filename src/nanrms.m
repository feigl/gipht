function nan_rms = nanrms(y)
% NAN_RMS calculate Root Mean Square of finite values of a vector or array
% 20180605 Kurt Feigl
%    iok = find(isfinite(y) == 1);
%    nan_rms = rms(y(iok));  
   nan_rms = rms(y(isfinite(y)));
end

