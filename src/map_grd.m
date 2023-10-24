function h = map_grd(varargin)
%function h = map_grd(grdfilename,cmap,SYM)
%function h = map_grd(grdfilename)
% map a GMT grid file named grdfilename and return graphics handle
% 20160814 Kurt Feigl
% 20170830 add SYM argument for plotting
% 20200129 TODO return INFO as optional output

%% parse input arguments
if nargin >= 1
    grdfilename = varargin{1};
else
    warning('missing arguments');
    h = nan;
    return;
end

% second argument is color table
if nargin >= 2
    cmap        = varargin{2};
else
    if contains(grdfilename,'amp') == true
        cmap = gray;
    else
        cmap = jet;
    end
end

% third argument is a structure
if nargin >= 3
    if isstruct(varargin{3}) == 1
        plot_symbols = 1;
        SYMS = varargin{3};
    else
        warning('Input argument SYM must be a structure');
        plot_symbols = 0;
    end
else
    plot_symbols = 0;
end

if nargin >= 4
    draw_grid = varargin{4}
else
    draw_grid = 0;
end


%% read input GMT grid file and its metadata
[xe,ye,IMAGE] = grdread4(grdfilename);
INFO = grdinfo4(grdfilename);

%% if coordinates are in UTM meters, then plot in kilometers
if isfield(INFO,'xname')
    if    contains(INFO.xname,'meters') == 1 ...
            || strcmp (INFO.yname,'m') == 1 ...
            || strcmp (INFO.yname,'m') == 1
        lengthfact = 1.e3; % plot in km
        %     xlab = strrep(INFO.xname,'in meters','[km]');
        %     ylab = strrep(INFO.yname,'in meters','[km]');
        xlab = 'Easting [km]';
        ylab = 'Northing [km]';
    else
        lengthfact = 1; % no scaling
        xlab = INFO.xname;
        ylab = INFO.yname;
    end
else
    lengthfact = 1; % no scaling
    xlab = 'Longitude [deg]';
    ylab = 'Latitude [deg]';
end
    

%% set up figure
figure; hold on;

% %% set up a symmetric colortable
% cmaxabs = nanmax(abs(colvec(IMAGE)));
% clim = [-cmaxabs,+cmaxabs];

%% set up a colortable
%clim = [nanmin(colvec(IMAGE)),nanmax(colvec(IMAGE))];
clim = [quantile(colvec(IMAGE),0.05),quantile(colvec(IMAGE),0.95)];
if abs(clim(1)-clim(2)) <= eps
    warning('mininum and maximum values are equal');
    %clim(1) = clim(1) - 0.1*abs(clim(1));
    %clim(2) = clim(2) + 0.1*abs(clim(2));
    clim(1) = -Inf;
    clim(2) = +Inf;
end

%% draw the image
%imagesc(xe/lengthfact,ye/lengthfact,IMAGE,clim);
%pcolor(xe/lengthfact,ye/lengthfact,IMAGE);
% set the transparency for values with NaN
ALPHA=ones(size(IMAGE));
ALPHA(isnan(IMAGE))=0.1;
imagesc(xe/lengthfact,ye/lengthfact,IMAGE,'AlphaData',ALPHA);

colormap(cmap);
axis xy;
axis equal;
axis tight;

if draw_grid == 1
    grid on;
    %Grid-line transparency, specified as a value in the range [0,1]. A
    %value of 1 means opaque and a value of 0 means completely transparent.
    set(gca,'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5],'GridAlpha',1,'Layer','top');
end

%% plot symbols if requested
if plot_symbols == 1
    for i=1:numel(SYMS.x)
        plot(SYMS.x(i)/lengthfact,SYMS.y(i)/lengthfact,SYMS.sym{i});
    end
end

%% add labels and title
xlabel(xlab);
ylabel(ylab);
tstr = sprintf('%s [%s]\n%s',INFO.title,grdfilename,INFO.description);
title(tstr,'Interpreter','None');

%% make color bar
c = colorbar;
if isfield(INFO,'zname')
    c.Label.String = INFO.zname;
end
h = gcf;

return

