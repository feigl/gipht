function r=rwrapm(x)
% wrap a value in radians onto [-pi,pi]
% NOTE: 2012-OCT-08
% >> rwrapm(pi)
% ans =
%    -3.1416
% >> rwrapm(-pi)
% ans =
%     3.1416
%r=2.0*pi*(x/pi/2.0 - round(x/pi/2.0));

% NOTE: 2012-OCT-08
r = angle(complex(cos(x),sin(x)));
% >> rwrapm(pi)
% ans =
%     3.1416
% >> rwrapm(-pi)
% ans =
%    -3.1416
return