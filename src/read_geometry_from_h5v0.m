close all
clear variables
%cd /System/Volumes/Data/mnt/t31/insar/SANEM/SENTINEL/T144d/mintpyRAMP/work
fname = 'inputs/geometryRadar.h5'
I=h5info(fname)
I.Datasets
nDataSets=numel(I.Datasets)
height=h5read(fname,'/height');
azimuthAngle=h5read(fname,'/azimuthAngle');
incidencAngle=h5read(fname,'/incidenceAngle');
latitude=h5read(fname,'/latitude');
longitude=h5read(fname,'/longitude');
shadowMask=h5read(fname,'/shadowMask');
slantRangeDistance=h5read(fname,'/slantRangeDistance');

return
for i=1:nDataSets
    name1=I.Datasets(i).Name
    G.(name1)=h5read(fname,sprintf('//%s',name1));
end

return

figure;
title(fname)
imagesc(Height);
axis xy
xlabel('Longitude');
ylabel('Latitude');
colormap(jet);
colorbar;
figure;
% E=h5read(fname,'/easting');
% N=h5read(fname,'/northing');
% plot(E/1.e3,N/1.e3,'.');
% xlabel('UTM Easting [km]');
% ylabel('UTM Northing [km]');
% title(fname)