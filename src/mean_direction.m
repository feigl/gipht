function thetabar = mean_direction(theta1)
% mean_direction(theta)
% given angles theta in radians,
% return mean direction in radians
% Kurt Feigl
%
% from Mardia 1972
% from Mardia and Jupp page 15
%
% This routine prunes NaN values
%
% 2011-JUN-22 update for third case
% according to N.I. Fisher equation 2.9, page 31
% Also handle special cases of Sbar == 0 || Cbar == 0
% 
% Examples:
%
% mean_direction([359,1,3]*pi/180)*180/pi
%
% ans =
%
%     1.0000
%
% mean_direction([360,0,0]*pi/180)*180/pi
% 
% ans =
% 
%      0


% prune NaN
iok = find(isfinite(theta1)==1);
theta = theta1(iok);

n = numel(theta);
if n > 0
%     Cbar = sum(cos(theta))/n;
%     Sbar = sum(sin(theta))/n;
    Cbar = sum(cos(theta));
    Sbar = sum(sin(theta));
else
    %warning('Number of data points is less than or equal to zero.');
    thetabar = NaN;
    return;
end
%xbar = angle(complex(Cbar/Rbar,Sbar/Rbar));
% Rbar = sqrt(Cbar^2 + Sbar^2);
% if Rbar <= 0.0 | Rbar > 1.0
%     thetabar = NaN;
% else
%     if Cbar >= 0
%         thetabar = atan(Sbar/Cbar);
%     else
%         thetabar = atan(Sbar/Cbar)+pi;
%     end
%     2011-JUN-22
% end
% Equation 2.9 does not define special case
if abs(Cbar) < 1.0e-15 || abs(Sbar) < 1.0e-15
    thetabar = 0.0;
elseif Sbar > 0 && Cbar > 0
    thetabar = atan(Sbar/Cbar);
elseif Cbar < 0
    thetabar = atan(Sbar/Cbar) + pi;
elseif Sbar < 0 && Cbar > 0
    thetabar = atan(Sbar/Cbar) + 2.0*pi;
else
    warning(sprintf('Undefined case Sbar = %10.4f Cbar = %10.4f\n',Sbar,Cbar));
    thetabar = NaN;
end

if isreal(thetabar) ~= 1 || isfinite(thetabar) ~= 1 
    warning(sprintf('Undefined case Sbar = %10.4f Cbar = %10.4f thetabar = %g\n',Sbar,Cbar, thetabar));
    thetabar = NaN;
end
    
return

