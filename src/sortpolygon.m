function [xout,yout] = sortpolygon(xin,yin)
%function [xout,yout] = sortpolygon(xin,yin)
% sort coordinates to make a closed polygon suitable for drawing with line function 
% Kurt Feigl 20140610
k = convhull(xin,yin);
xout = xin(k);
yout = yin(k);
xout(end+1) = xout(1);
yout(end+1) = yout(1);
end

