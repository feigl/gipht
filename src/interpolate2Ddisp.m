function [upts,vpts,wpts] = interpolate2Ddisp(PST,xobs,yobs,rnode,znode,hnode,wnode)



% establish polar coordinates
% angle CCW from East
theta = atan2(yobs - PST.ycen,xobs - PST.xcen);

% radial distance from center
rdist = sqrt((xobs - PST.xcen).^2 + (yobs - PST.ycen).^2);

% evaluate functions at ptsput locations
% dhpts  = interp1(rnode,dhnode,rdist,'nearest','extrap'); % horizontal component
% dwpts  = interp1(rnode,dwnode,rdist,'nearest','extrap'); % vertical component
%
hpts  = interp1(rnode,hnode,rdist,'cubic','extrap'); % horizontal component
wpts  = interp1(rnode,wnode,rdist,'cubic','extrap'); % vertical component

% ibad = find(isfinite(dhpts) ~= 1);
% if numel(ibad) > 0
%     warning('Found %d NaN values\n',numel(ibad));
%     dhpts(ibad) = 0;
% end
% ibad = find(isfinite(dwpts) ~= 1);
% if numel(ibad) > 0
%     warning('Found %d NaN values\n',numel(ibad));
%     wpts(ibad) = 0;
% end

% incremental displacement between t1 and t2
upts  = hpts .* cos(theta);       % east component
vpts  = hpts .* sin(theta);       % north component



end

