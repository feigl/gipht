function GEOM = read_geometry_from_h5(fname,doplots)
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

if exist('doplots','var') ~= true
    doplots=true;
end
for i=1:nDataSets
    name1=I.Datasets(i).Name
    
    D = h5read(fname,sprintf('//%s',name1)); 
    if strcmpi(name1,'date')
        yyyymmdd=D;
    end
    if isnumeric(D)
        nd=ndims(D)
        switch nd
            case 1
                GEOM.(name1)=D;
            case 2
                GEOM.(name1)=flipud(transpose(D));
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
            case 3
                [nx,ny,nt] = size(D)
                for it=1:nt
                    if exist('yyyymmdd','var')
                        slice1=sprintf('T%8s',yyyymmdd(it));
                    else
                        slice1=sprintf('slice%04d',it);
                    end
                    GEOM.(slice1)=flipud(transpose(squeeze(D(:,:,it))));
                    if doplots
                        figure;
                        imagesc(GEOM.(slice1));
                        axis xy
                        xlabel('X index');
                        ylabel('Y index');
                        colormap(jet);
                        title(sprintf('%s\n%s',fname,name1),'Interpreter','None');
                        colorbar;
                    end
                end
            otherwise
                warning(sprintf('unknown number of dimensions nd = %d\n',nd));
        end
    else
        GEOM.(name1)=D;
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


