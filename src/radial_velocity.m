function V = radial_velocity(params,xobs,yobs)
% given parameters, calculate 3 components of velocity
% input parameters param:
%  xcen == Rad Vel Center Easting in m
%  ycen == Rad Vel Center Northing in m
%  ucen == Rad Vel Evel in m per yr
%  vcen == Rad Vel Nvel in m per yr
%  wcen == Rad Vel Uvel in m per yr
%  dvdr == Rad Vel dVdR in inverse yr
%  xobs == easting coordinate of observation point in meters
%  yobs == northing coordinate of observation point in meters
% output
% [ve, vn, vu] == [eastward,northward,upward] components of velocity in meters/year

% check parameters
if numel(params) == 6
    xcen=params(1);
    ycen=params(2);
    ucen=params(3);
    vcen=params(4);
    wcen=params(5);
    dvdr=params(6);
else
    error(sprintf('Number of parameters %d is not equal to %d\n',numel(params),6));   
end

% check data
ndata = numel(xobs);
if numel(yobs) == ndata
     V = zeros(3,ndata);
else
    error('dimension mismatch');   
end

% radial distance
rdist = rowvec(sqrt((xobs-xcen).^2 + (yobs-ycen).^2));
% fprintf(1,'Extrema of rdist %e %e\n',min(rdist),max(rdist));
%rdist = reshape(rdist,numel(rdist),1);
%rmax  = nanmax(rdist);

% velocity decays from max at center
V(1,:) = ucen*rdist*dvdr;
V(2,:) = vcen*rdist*dvdr;
V(3,:) = wcen*rdist*dvdr;
% fprintf(1,'Extrema of V east  %e %e\n',min(V(1,:)),max(V(1,:)));
% fprintf(1,'Extrema of V north %e %e\n',min(V(2,:)),max(V(2,:)));
% fprintf(1,'Extrema of V up    %e %e\n',min(V(3,:)),max(V(3,:)));
% 
% figure
% %quiver(colvec(xobs),colvec(yobs),colvec(V(1,:)),colvec(V(3,:)));
% plot(colvec(xobs),colvec(V(3,:)),'k+');
% title(mfilename);

return

