function [U1,U2,U3]=fcn_yangM(as,x,y,matrl,tp)
% Calculate range displacement 
mu=matrl(2);
nu=matrl(3);

% Store some commonly used parameters
coeffs(1)=1/(16*mu*(1-nu));
coeffs(2)=3-4*nu;
coeffs(3)=4*(1-nu)*(1-2*nu);

[ix,iy]=size(x);
U1r=zeros(ix,iy);
U2r=zeros(ix,iy);
U1=zeros(ix,iy);
U2=zeros(ix,iy);
U3=zeros(ix,iy);

% explicitly assign source parameters
%as
 xs    = as(1);     % center x
 ys    = as(2);     % center y
 z0    = as(3);     % center depth (positive)
 P     = as(4);     % excess pressure, mu*10^(-5) Pa
 a     = as(5);     % major axis, km
 b     = as(6);     % minor axis, km
 
 % 20130419 catch case where spheroid is oblate, rather than spheroid
 if b > a 
     warning(sprintf('Semiminor axis b (%f) is greater than semimajor axis a (%f)\n',b,a));
     return;
 end
 phi   = as(7)*pi/180.0;     % strike, rad  (0-2*pi)
 theta = as(8)*pi/180.0;     % plunge, rad  (0-pi)
 xn=x-xs;
 yn=y-ys;
 e_theta(1)=sin(theta);
 e_theta(2)=cos(theta);
 cosp=cos(phi);
 sinp=sin(phi);
 c=sqrt(a^2-b^2);
 minx=min(min(x));
 maxx=max(max(x));
 miny=min(min(y));
 maxy=max(max(y));
% if xs < minx | xs > maxx | ys < miny | ys > maxy % source is outside the grid
%  P=0;
% end

% Speroid quantities
 [sph]=spheroid(a,b,c,matrl,phi,theta,P);

% Rotate points
 xp=xn*cosp + yn*sinp;
 yp=yn*cosp - xn*sinp;

% Calculate model at integration limits
 xi=c;
 [Up1,Up2,Up3]=yang(sph,xi,z0,xp,yp,0,matrl,e_theta,coeffs,tp);
 xi=-xi;
 [Um1,Um2,Um3]=yang(sph,xi,z0,xp,yp,0,matrl,e_theta,coeffs,tp);

% Sum
 U1r=-Up1+Um1;
 U2r=-Up2+Um2;
% Rotate horiz. displacements back to the orig. coordinate system:
 U1=U1r*cosp-U2r*sinp+U1;
 U2=U1r*sinp+U2r*cosp+U2;
 U3=Up3-Um3+U3;
 
%  ibad = find(isfinite(U1)==0);
%  if numel(ibad) > 0 
%      warning ('found Nan numbers in U1');
%      U1(ibad) = 0;U2(ibad) = 0;U3(ibad) = 0;
%  end
%  ibad = find(isfinite(U2)==0);
%  if numel(ibad) > 0 
%      warning ('found Nan numbers in U2');
%      U1(ibad) = 0;U2(ibad) = 0;U3(ibad) = 0;
%  end
% 
%  ibad = find(isfinite(U3)==0);
%  if numel(ibad) > 0 
%      warning ('found Nan numbers in U3');
%      U1(ibad) = 0;U2(ibad) = 0;U3(ibad) = 0;
%  end
 



