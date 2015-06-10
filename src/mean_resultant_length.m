function r = mean_resultant_length (theta1)
% mean_resultant_length (theta1)
% given angles theta in radians,
% return resltant_length Rbar in dimesionless units
% Kurt Feigl
% from Mardia 1972
%

% prune NaN
iok = find(isfinite(theta1)==1);
theta=theta1(iok);

n=numel(theta);
if n > 0
    cbar=sum(cos(theta))/n;
    sbar=sum(sin(theta))/n;
    r=sqrt(cbar^2 + sbar^2);
    
    % sanity check
    if r < 0
        r = 0;
    end
    if r > 1
        r = 1;
    end
else
    r = NaN;
end

return
