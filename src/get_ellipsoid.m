function [x,y,z] = get_ellipsoid(xc,yc,zc,a,b,c,azimuth,plunge)
%function [x,y,z] = get_ellipsoid(xc,yc,zc,a,b,c,azimuth,plunge)
% calculate locus of points on an elllipsoid with three axes
%
% inputs
%     
% (xc,yc,zc) center in meters (a is semi-major, b is semi-minor)
% azimuth in degrees clockwise from north
% plunge in degrees downward from horizontal

npoints = 1000;

% convert from degrees to radians
azimuth = azimuth * pi/180.0;
plunge  = plunge  * pi/180.0;

[x,y,z]=ellipsoid(xc,yc,zc,a,b,c,npoints);

% figure;axis(1.1*[-a, a, -a, a, -a, a]);hold on;
% xlabel('easting (m)');ylabel('northing (m)');zlabel('upward (m)');
% surf(x,y,z)
x=colvec(x);
y=colvec(y);
z=colvec(z);
npoints2 = numel(x);

% rotate about vertical axis by complement of azimuth
%azimuth = 75 * pi /180.;
%azimuth = 0;
% azimuth = 90 * pi /180.;

R=rotationmat3D(pi/2 - azimuth,[0,0,1]);
for i=1:npoints2
    X(1) = x(i)-xc;
    X(2) = y(i)-yc; 
    X(3) = z(i)-zc;
    X2 = R * X';
    x(i) = X2(1)+xc;
    y(i) = X2(2)+yc;
    z(i) = X2(3)+zc;
end
   
% figure;axis(1.1*[-a, a, -a, a, -a, a]);hold on;
% xlabel('easting (m)');ylabel('northing (m)');zlabel('upward (m)');
% surf(reshape(x,npoints+1,npoints+1) ...
%     ,reshape(y,npoints+1,npoints+1) ...
%     ,reshape(z,npoints+1,npoints+1));

% rotate about horizontal axis by plunge
%plunge = 90 * pi /180.; 
%plunge = 0;
R=rotationmat3D(plunge,[sin(pi/2-azimuth),cos(pi/2-azimuth),0]);
for i=1:npoints2
    X(1) = x(i)-xc;
    X(2) = y(i)-yc; 
    X(3) = z(i)-zc;
    X2 = R * X';
    x(i) = X2(1)+xc;
    y(i) = X2(2)+yc;
    z(i) = X2(3)+zc;
end

% figure;axis(1.1*[-a, a, -a, a, -a, a]);hold on;
% xlabel('easting (m)');ylabel('northing (m)');zlabel('upward (m)');
% surf(reshape(x,npoints+1,npoints+1) ...
%     ,reshape(y,npoints+1,npoints+1) ...
%     ,reshape(z,npoints+1,npoints+1));
% 
% % draw slice at center depth
% figure;axis(1.1*[-a, a, -a, a]);hold on
% xlabel('easting (m)');ylabel('northing (m)');
% islice = find(abs(z-zc)<5.);
% plot(x(islice),y(islice),'r.');
% xmin=min(x(islice));
% xmax=max(x(islice));
% ymin=min(y(islice));
% ymax=max(y(islice));
% xrect = [xmin xmax xmax xmin xmin];
% yrect = [ymin ymin ymax ymax ymin];
% plot(xrect,yrect,'k-');

return




   
