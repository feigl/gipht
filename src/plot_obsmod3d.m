function h = plot_obsmod3d(DST)
%function h = plot_obscal3(DST)       
% plot the observed and modeled values in 3 dimensions
figure; 

% interpolate and plot exact modeled values
[XI,YI] = meshgrid(linspace(min(DST.x),max(DST.x)),linspace(min(DST.y),max(DST.y)));
ZI = griddata2(DST.x,DST.y,DST.phamod,XI,YI);
surf(XI/1.e3,YI/1.e3,ZI);
hold on;
title('observed values (black dots) and modeled values (colored surface)');
xlabel('Easting [km]');
ylabel('Northing [km]');
zlabel('Range [radians]');

% plot observed values
%plot3(DST.x/1e3,DST.y/1e3,DST.phaobs,'ko','MarkerFaceColor','k');
for i = 1:numel(DST.phaobs)
    plot3([DST.x(i); DST.x(i)]/1e3 , [DST.y(i); DST.y(i)]/1e3,[0; DST.phaobs(i)],'k-');
    plot3(DST.x(i)/1e3,DST.y(i)/1e3,DST.phaobs(i),'ko','MarkerFaceColor','k');
end

axis tight;

h = gcf;

return
end

