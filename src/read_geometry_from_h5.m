function GEOM = read_geometry_from_h5(fname,varargin)
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
% isgeographic == true for geographic
% 2023/01/01 Kurt Feigl
% 2024/04/29 Kurt Feigl
% 2024/09/09

nf=0;
% set default values
if nargin < 2
    doplots=true;
elseif nargin < 3
    doplots=varargin{1};
    isGeographic=false;
else
    doplots=varargin{1};
    isGeographic=varargin{2};
end
doplots
isGeographic

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
            fprintf(1,'%s (%s) \t = %.10g\n',name1,name2,value1);
            GEOM.atts.(name2)=value1;
        end
    end
    % make vectors of coordinates
    %     GEOM.x_vec=[GEOM.atts.X_FIRST - GEOM.atts.X_STEP/2. : GEOM.atts.X_STEP : GEOM.atts.X_FIRST + (GEOM.atts.WIDTH-1)  * GEOM.atts.X_STEP + GEOM.atts.X_STEP/2.];
    %     GEOM.y_vec=[GEOM.atts.Y_FIRST - GEOM.atts.Y_STEP/2. : GEOM.atts.Y_STEP : GEOM.atts.Y_FIRST + (GEOM.atts.LENGTH-1) * GEOM.atts.Y_STEP + GEOM.atts.Y_STEP/2.];
    %GEOM.x_vec=[GEOM.atts.X_FIRST  : GEOM.atts.X_STEP : GEOM.atts.X_FIRST + (GEOM.atts.WIDTH-1)  * GEOM.atts.X_STEP ] + GEOM.atts.X_STEP/2. ;
    %GEOM.y_vec=[GEOM.atts.Y_FIRST  : GEOM.atts.Y_STEP : GEOM.atts.Y_FIRST + (GEOM.atts.LENGTH-1) * GEOM.atts.Y_STEP ] + GEOM.atts.Y_STEP/2. ;
    %     GEOM.x_vec=GEOM.atts.X_FIRST + [0 : GEOM.atts.X_STEP : (GEOM.atts.WIDTH-1)   * GEOM.atts.X_STEP ] + GEOM.atts.X_STEP/2. ;
    %     GEOM.y_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] + GEOM.atts.Y_STEP/2. ;
    %     GEOM.x_vec=GEOM.atts.X_FIRST + [0 : GEOM.atts.X_STEP : (GEOM.atts.WIDTH-1)   * GEOM.atts.X_STEP ] ;
    %     GEOM.y_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] ;

    GEOM.x_vec=GEOM.atts.X_FIRST + [0 : GEOM.atts.X_STEP : (GEOM.atts.WIDTH-1)   * GEOM.atts.X_STEP ] ;
    % 20240814
    %     if GEOM.atts.Y_STEP > 0
    %         GEOM.y_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] ;
    %     else
    %         GEOM.y_vec=GEOM.atts.Y_FIRST + [(GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP : GEOM.atts.Y_STEP : 0 ] ;
    %     end
    % 20240814 applies even if step size is negative
    GEOM.y_vec=GEOM.atts.Y_FIRST + [0 : GEOM.atts.Y_STEP : (GEOM.atts.LENGTH-1)  * GEOM.atts.Y_STEP ] ;

    % make a grid and re-arrange from Python convention to Matlab convention
    [GEOM.x,GEOM.y] = meshgrid(GEOM.x_vec,GEOM.y_vec);
    % GEOM.x = flipud(transpose(GEOM.x));
    % GEOM.y = flipud(transpose(GEOM.y));
    GEOM.x = flipud(GEOM.x);
    GEOM.y  = flipud(GEOM.y);
    % sort the coordinates
    GEOM.x_vec = sort(GEOM.x_vec);
    GEOM.y_vec = sort(GEOM.y_vec);
end
%GEOM

for i=1:nDataSets
    name1=I.Datasets(i).Name

    clear A;
    fprintf(1,'Reading data set %s\n',name1)
    try
        A = h5read(fname,sprintf('//%s',name1));
    catch ME
        ME.message
        for ii=1:numel(ME.stack)
            ME.stack(ii)
        end
        warning('failed');
        A=nan;
    end
    if exist('A','var')==1
        if isnumeric(A)
            A = double(A);
            ndim=ndims(A)
            if ndim == 2
                [nx,ny] = size(A)
                if nx > 1 && ny > 1
                    GEOM.(name1)=flipud(transpose(A));
                else
                    warning('strange dimensions')
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

            if ~isfield(GEOM,'lon_vec') % || ~isfield(GEOM,'x')
                [nlon,mlon] = size(GEOM.x)
                %mlon=nx
                GEOM.x_vec=linspace(nanmin(colvec(GEOM.x)),nanmax(colvec(GEOM.x)),mlon);
            end
            if ~isfield(GEOM,'lat_vec') % || ~isfield(GEOM,'y')
                [nlat,mlat] = size(GEOM.y)
                %mlat=ny
                GEOM.y_vec=linspace(nanmin(colvec(GEOM.y)),nanmax(colvec(GEOM.y)),mlat);
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
                                imagesc(GEOM.x_vec,GEOM.y_vec,squeeze(GEOM.(name1)(:,:,ii)));
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
                            imagesc(GEOM.x_vec,GEOM.y_vec,GEOM.(name1));
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
            % rename fields
            if isGeographic && isfield(GEOM,'x') && isfield(GEOM,'y')
                GEOM.longitude=GEOM.x;    GEOM=rmfield(GEOM,'x');
                GEOM.latitude =GEOM.y;    GEOM=rmfield(GEOM,'y');
                GEOM.lon_vec  =GEOM.x_vec;GEOM=rmfield(GEOM,'x_vec');
                GEOM.lat_vec  =GEOM.y_vec;GEOM=rmfield(GEOM,'y_vec');
            end
         else
             % A is not numeric
             if strcmp(name1,'date')
                 % A is an array of strings
                 % size(A)
                 ndates=length(A)
                 GEOM.date=A;
                 for ii=1:ndates
                     %Astring=sprintf('%8d\n',A(ii));
                     Astring=A(ii);
                     % convert date, ignoring time of day
                     dt=datetime(Astring,'InputFormat','yyyyMMdd','TimeZone','UTC','format',['yyyy-MM-dd''T''HH:mmXXXZZZZ']);
                     GEOM.Time(ii)=dt;
                 end
             else
                 warning('data set %s is not numeric',name1)
                 GEOM.(name1)=A;
             end
        end
    else
        GEOM.(name1)=nan;
    end
end
return
end

