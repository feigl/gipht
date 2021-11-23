function GEOM = read_geometry_from_h5(fname)
% read an H5 file from MINTPY and return structure with arrays containing geometry
%cd /System/Volumes/Data/mnt/t31/insar/SANEM/SENTINEL/T144d/mintpyRAMP/work
%fname = 'inputs/geometryRadar.h5'
I=h5info(fname)
I.Datasets
nDataSets=numel(I.Datasets)
% height=h5read(fname,'/height');
% azimuthAngle=h5read(fname,'/azimuthAngle');
% incidencAngle=h5read(fname,'/incidenceAngle');
% latitude=h5read(fname,'/latitude');
% longitude=h5read(fname,'/longitude');
% shadowMask=h5read(fname,'/shadowMask');
% slantRangeDistance=h5read(fname,'/slantRangeDistance');

doplots=true;
for i=1:nDataSets
    name1=I.Datasets(i).Name
    GEOM.(name1)=flipud(transpose(h5read(fname,sprintf('//%s',name1))));
    
    if doplots
        figure;        
        imagesc(GEOM.(name1));
        axis xy
        xlabel('X index');
        ylabel('Y index');
        colormap(jet);
        title(sprintf('%s\n%s',fname,name1),'Interpreter','None');
        colorbar;
    end    
end
if isfield(GEOM,'longitude')
   [nlon,mlon] = size(GEOM.longitude);
   GEOM.lon_vec=linspace(nanmin(colvec(GEOM.longitude)),nanmax(colvec(GEOM.longitude)),mlon);
end
if isfield(GEOM,'latitude')
   [nlat,mlat] = size(GEOM.latitude);
   GEOM.lat_vec=linspace(nanmin(colvec(GEOM.latitude)),nanmax(colvec(GEOM.latitude)),mlat);
end
return
end


