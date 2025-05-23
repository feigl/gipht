<<<<<<< HEAD
function GEOM = read_geometry_from_h5(fname,doplots)
=======
function GEOM = read_geometry_from_h5(fname,varargin)
>>>>>>> c0a77c90e8512fa2eda6314b5513fb71fb565958
% read an H5 file from MINTPY and return structure with arrays containing geometry
% sources of information
% https://github.com/insarlab/MintPy-tutorial/blob/main/smallbaselineApp_aria.ipynb
% https://mintpy.readthedocs.io/en/latest/FAQs/
%
% % Example for WHOLESCALE
%cd /System/Volumes/Data/mnt/t31/insar/SANEM/SENTINEL/T144d/mintpyRAMP/work
%fname = 'inputs/geometryRadar.h5'
% display
%h5disp(fname)
% height=h5read(fname,'/height');
% azimuthAngle=h5read(fname,'/azimuthAngle');
% incidencAngle=h5read(fname,'/incidenceAngle');
% latitude=h5read(fname,'/latitude');
% longitude=h5read(fname,'/longitude');
% shadowMask=h5read(fname,'/shadowMask');
% slantRangeDistance=h5read(fname,'/slantRangeDistance');
% 2023/01/01 Kurt Feigl

<<<<<<< HEAD
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
=======
nf=0;
% set default values
if nargin < 2
    doplots=true;
else
    doplots=varargin{1};
end

%% get attributes from ARIA file
I=h5info(fname)
I.Datasets
I.Datasets.Name
I.Attributes
nDataSets=numel(I.Datasets)

if isfield(I,'Attributes')
    GEOM.Attributes=I.Attributes;
    nAttributes = length(GEOM.Attributes);
    % convert numeric attributes from text to scalar values 
    % this is clumsy and slow, but I don't know how to do better
    [names{1:nAttributes,1}]=deal(GEOM.Attributes.Name);
    [values{1:nAttributes,1}]=deal(GEOM.Attributes.Value);
    for ii=1:nAttributes
        value1=str2double(values{ii});
        name1 = names{ii};
        name2 = matlab.lang.makeValidName(name1);
        if isfinite(value1)
            %fprintf(1,'%s (%s) \t = %.10g\n',name1,name2,value1);
            GEOM.atts.(name2)=value1;
        end
    end
    % make vectors of coordinates
%     GEOM.lon_vec=[GEOM.atts.X_FIRST - GEOM.atts.X_STEP/2. : GEOM.atts.X_STEP : GEOM.atts.X_FIRST + (GEOM.atts.WIDTH-1)  * GEOM.atts.X_STEP + GEOM.atts.X_STEP/2.];
%     GEOM.lat_vec=[GEOM.atts.Y_FIRST - GEOM.atts.Y_STEP/2. : GEOM.atts.Y_STEP : GEOM.atts.Y_FIRST + (GEOM.atts.LENGTH-1) * GEOM.atts.Y_STEP + GEOM.atts.Y_STEP/2.];
    %GEOM.lon_vec=[GEOM.atts.X_FIRST  : GEOM.atts.X_STEP : GEOM.atts.X_FIRST + (GEOM.atts.WIDTH-1)  * GEOM.atts.X_STEP ] + GEOM.atts.X_STEP/2. ;
    %GEOM.lat_vec=[GEOM.atts.Y_FIRST  : GEOM.atts.Y_STEP : GEOM.atts.Y_FIRST + (GEOM.atts.LENGTH-1) * GEOM.atts.Y_STEP ] + GEOM.atts.Y_STEP/2. ;
%     GEOM.lon_vec=GEOM.atts.X_FIRST + [0 : GEOM.atts.X_STEP : (GEOM.atts.WIDTH-1)   * GEOM.atts.X_STEP ] + GEOM.atts.X_STEP/2. ;
%     GEOM.lat_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] + GEOM.atts.Y_STEP/2. ;
    GEOM.lon_vec=GEOM.atts.X_FIRST + [0 : GEOM.atts.X_STEP : (GEOM.atts.WIDTH-1)   * GEOM.atts.X_STEP ] ;
    GEOM.lat_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] ;

    % make a grid and re-arrange from Python convention to Matlab convention
    [GEOM.longitude,GEOM.latitude] = meshgrid(GEOM.lon_vec,GEOM.lat_vec);
    % GEOM.longitude = flipud(transpose(GEOM.longitude));
    % GEOM.latitude = flipud(transpose(GEOM.latitude));
    GEOM.longitude = flipud(GEOM.longitude);
    GEOM.latitude  = flipud(GEOM.latitude);
    % sort the coordinates
    GEOM.lon_vec = sort(GEOM.lon_vec);
    GEOM.lat_vec = sort(GEOM.lat_vec);
end

for i=1:nDataSets
    name1=I.Datasets(i).Name
    clear A;
    fprintf(1,'Reading data set %s\n',name1)
    A = h5read(fname,sprintf('//%s',name1));
    if ~isnumeric(A)
        A = str2double(A);
    end
    ndim=ndims(A)
    if ndim == 2
        [nx,ny] = size(A)
        if nx > 1 && ny > 1
            GEOM.(name1)=flipud(transpose(A));            
        else
            size(A)
            if strcmp(name1,'date')
                ndates=length(A);
                GEOM.date=A;
                for ii=1:ndates
                    Astring=sprintf('%8d\n',A(ii));
                    % convert date, ignoring time of day
                    dt=datetime(Astring,'InputFormat','yyyyMMdd','TimeZone','UTC','format',['yyyy-MM-dd''T''HH:mmXXXZZZZ']);
                    GEOM.Time(ii)=dt;
                end
            else
                GEOM.(name1)=A;
            end
        end
    elseif ndim == 3 && strcmp(name1,'timeseries')
        [ny,nx,nt] = size(A);
        for ii=1:nt
            GEOM.(name1)(:,:,ii)=flipud(transpose(squeeze(A(:,:,ii))));
        end
    else
        size(A)
        GEOM.(name1)=A;
    end

    if ~isfield(GEOM,'lon_vec')
        [nlon,mlon] = size(GEOM.longitude);
        GEOM.lon_vec=linspace(nanmin(colvec(GEOM.longitude)),nanmax(colvec(GEOM.longitude)),mlon);
    end
    if ~isfield(GEOM,'lat_vec')
        [nlat,mlat] = size(GEOM.latitude);
        GEOM.lat_vec=linspace(nanmin(colvec(GEOM.latitude)),nanmax(colvec(GEOM.latitude)),mlat);
    end

    if doplots    
        switch name1
            case 'date'
                nf=nf+1;figure;
                bar(GEOM.(name1),ones(size(GEOM.(name1))));
                xlabel('date');
            case 'timeseries'
                fprintf(1,'plotting %d time series\n',nt);
                for ii=1:nt
                     nf=nf+1;figure;
                    if isfield(GEOM,'lon_vec') && isfield(GEOM,'lat_vec')
                        imagesc(GEOM.lon_vec,GEOM.lat_vec,squeeze(GEOM.(name1)(:,:,ii)));
                        xlabel('longitude');
                        ylabel('latitude');
                    else
                        imagesc(GEOM.(name1)(:,:,ii));
                        xlabel('X index');
                        ylabel('Y index');
                    end
                    axis xy; axis image;
                    colormap(jet);
                    colorbar;
                    title(sprintf('%s\n%s epoch %d date %s',fname,name1,ii,'date'),'Interpreter','None');
                end
            case {'longitude','latitude','azimuthAngle','height','incidenceAngle','slantRangeDistance','waterMask','velocity','velocityStd'}
                 nf=nf+1;figure;
                if isfield(GEOM,'lon_vec') && isfield(GEOM,'lat_vec')
                    imagesc(GEOM.lon_vec,GEOM.lat_vec,GEOM.(name1));
                    xlabel('longitude');
                    ylabel('latitude');                   
                else
                    imagesc(GEOM.(name1));
                    xlabel('X index');
                    ylabel('Y index');
                end
                axis xy; axis image;
                colormap(jet);
                colorbar;
                title(sprintf('%s\n%s',fname,name1),'Interpreter','None');
            otherwise
                 nf=nf+1;figure;
                plot(GEOM.(name1),'ro-');
                axis xy;
                xlabel('X index');
                ylabel('value');
                title(sprintf('%s\n%s',fname,name1),'Interpreter','None');
        end
        exportgraphics(gcf,sprintf('%s_%03d.png',mfilename,nf),'Resolution',1200);
    end  
>>>>>>> c0a77c90e8512fa2eda6314b5513fb71fb565958
end
return
end


