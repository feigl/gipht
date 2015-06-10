function r=rwrapm2(x)
% wrap a value in radians onto [-pi,pi]
% for values very near +pi, return -pi, such that
% rwrapm2(-pi)
% ans =
%   -3.141592653589793
% rwrapm2(pi)
% ans =
%   -3.141592653589793
% Kurt Feigl 2014-01-06
r = angle(complex(cos(x),sin(x)));
ipi = find(pi - r < 2.0*pi/256.0);
r(ipi) = r(ipi) - 2.0 * pi;    
return
end
