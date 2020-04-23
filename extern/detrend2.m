function  [zo,r0,p] = detrend2(x,y,zi)
% DETREND2 Removing mean and slope in 2-d set.
%	[ZO,R0,P] = DETREND2(X,Y,ZI) calculates
%	and removes mean and slope of the dataset
%	with coordinates X, Y, Z.
%	Returns "detrended" value ZO so that
%	ZI = P(1)+(X-R0(1))*P(2)+(Y-R0(2))*P(3)+ZO.

%  Copyright (c) Kirill K. Pankratov
%	kirill@plume.mit.edu
%	05/20/95	

 % Handle input ......................
if nargin==0, help detrend2, return, end
if nargin<3
  error('  Not enough input arguments')
end

 % Make all coordinates column vectors
sz = size(zi);
x = x(:);
y = y(:);
zi = zi(:);

 % Find mean values
r0(1) = mean(x);
r0(2) = mean(y);
z0 = mean(zi);

 % Extract mean values
x = x-r0(1);
y = y-r0(2);
zo = zi-z0(1);

 % Find slope p = [zmean dx dy]
p = [x y]\zo;
p = [z0 p'];

 % Extract slope
zo = zo-x*p(2)-y*p(3);
zo = reshape(zo,sz(1),sz(2));
