function h = map_grd(grdfilename,cmap)
%function h = map_grd(grdfilename)
% map a GMT grid file named grdfilename and return graphics handle
% 20160814 Kurt Feigl
figure;
[xe,ye,IMAGE] = grdread3(grdfilename);
INFO = grdinfo3(grdfilename);

% if coordinates are in UTM meters, then plot in kilometers
if numel(strfind(INFO.xname,'meters')) > 0
    xlab = strrep(INFO.xname,'in meters','[km]');
    xp = xe/1.e3;
else
    xlab = INFO.xname;
    xp = xe;
end
if numel(strfind(INFO.yname,'meters')) > 0
    ylab = strrep(INFO.yname,'in meters','[km]');
    yp = ye/1.e3;
else
    ylab = INFO.yname;
    yp = ye;
end

% set up a symmetric colortable
cmaxabs = max(abs(colvec(IMAGE)));
clim = [-cmaxabs,+cmaxabs];

% draw the image
imagesc(xp,yp,IMAGE,clim);
colormap(cmap);
axis xy;
axis equal;
axis tight;
xlabel(xlab);
ylabel(ylab);
tstr = sprintf('%s\n[%s]',INFO.title,grdfilename);
title(tstr,'Interpreter','None');
c = colorbar;
c.Label.String = INFO.zname;
h = gcf;

return

