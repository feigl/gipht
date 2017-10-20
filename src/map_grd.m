function h = map_grd(varargin)
%function h = map_grd(grdfilename,cmap,SYM)
%function h = map_grd(grdfilename)
% map a GMT grid file named grdfilename and return graphics handle
% 20160814 Kurt Feigl
% 20170830 add SYM argument for plotting

%% parse input arguments
if nargin >= 1
    grdfilename = varargin{1};
else
    warning('missing arguments');
    h = nan;
    return;
end

if nargin >= 2
    cmap        = varargin{2};
else
    cmap = jet;
end

if nargin == 3
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

%% read input GMT grid file and its metadata
[xe,ye,IMAGE] = grdread3(grdfilename);
INFO = grdinfo3(grdfilename);

%% if coordinates are in UTM meters, then plot in kilometers
if contains(INFO.xname,'meters') == 1 || contains (INFO.yname,'meters') == 1
    lengthfact = 1.e3; % plot in km
    xlab = strrep(INFO.xname,'in meters','[km]');
    ylab = strrep(INFO.yname,'in meters','[km]');   
else
    lengthfact = 1; % no scaling
    xlab = INFO.xname;
    ylab = INFO.yname; 
end
    

%% set up figure
figure; hold on;

% %% set up a symmetric colortable
% cmaxabs = nanmax(abs(colvec(IMAGE)));
% clim = [-cmaxabs,+cmaxabs];

%% set up a colortable
clim = [nanmin(colvec(IMAGE)),nanmax(colvec(IMAGE))];

%% draw the image
imagesc(xe/lengthfact,ye/lengthfact,IMAGE,clim);
colormap(cmap);
axis xy;
axis equal;
axis tight;

%% plot symbols if requested
if plot_symbols == 1
    for i=1:numel(SYMS.x)
        plot(SYMS.x(i)/lengthfact,SYMS.y(i)/lengthfact,SYMS.sym{i});
    end
end

%% add labels and title
xlabel(xlab);
ylabel(ylab);
tstr = sprintf('%s\n[%s]',INFO.title,grdfilename);
title(tstr,'Interpreter','None');

%% make color bar
c = colorbar;
c.Label.String = INFO.zname;
h = gcf;

return

