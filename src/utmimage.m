function varargout = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax...
   ,titlestr,cornerlabel,climits,dotx,doty,ctab,mysym,marksize...
   ,drawxlabels,drawylabels,datelabel...
   ,idatatype,datalabel,cbar)
% function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax...
%    ,titlestr,cornerlabel,climits,dotx,doty,ctab,mysym,marksize...
%    ,drawxlabels,drawylabels,datelabel...
%    ,idatatype,datalabel)
% 2021/06/22 Kurt Feigl update to handle list of pixels
% 2021/09/02 Kurt Feigl 
%  if climit has more than 2 elements, use it to stretch color table

%% initialize 
% default is to return an error
istat = -1;
Vstretch=nan(size(pixarr));

if nargin < 8 || exist('climit','var') == 0 || exist('dotx','var') == 0 || exist('doty','var') == 0
   climit  = [-0.5 0.5];
   dotx = 0;
   doty = 0;
end
if nargin < 11 || exist('ctab','var') == 0
   ctab = 'jet';
end
if nargin < 12 || exist('mysym','var') == 0
   mysym = 'ko';  % default is draw a black circle, but do not connect
end
if nargin < 13 || exist('marksize','var') == 0
   marksize = 5;
end
if nargin < 14 || exist('drawxlabels','var') == 0
   drawxlabels=0;
end
if nargin < 15 || exist('drawylabels','var') == 0
   drawylabels=0;
end
if nargin < 16 || exist('datelabel','var') == 0
   datelabel='';
end
if nargin < 17 || exist('idatatype','var') == 0
   idatatype = 0; % wrapped phase is default
end
if nargin < 18 || exist('cbar','var') == 0
   cbar = 0; % default is to NOT draw a colorbar
end

%labelcolor = [0.5 0.5 0.5]; % gray
labelcolor = [0.0 0.0 0.0]; % black

%% make this routine robust
if isfinite(xutmmin) == 0
    xutmmin = 0;
end
if isfinite(yutmmin) == 0
    yutmmin = 0;
end
if isfinite(xutmmax) == 0
    xutmmax = 1;
end
if isfinite(yutmmax) == 0
    yutmmax = 1;
end

if xutmmax - xutmmin < 1.0
    xutmmax = xutmmin + 1;
end
if yutmmax - yutmmin < 1.0
    yutmmax = yutmmin + 1;
end

% if idatatype == 2
%     climit(1) = nanmin(colvec(pixarr));
%     climit(2) = nanmax(colvec(pixarr));
% else
%     % parse levels and extrema
%     iok = find(isfinite(climit)==1);
%     climit = climit(iok);
%     nlevels = numel(climit);
%     if nlevels < 2
%         cuts = NaN;
%         climit(1) = -0.5;
%         climit(2) = +0.5;
%     elseif nlevels == 2
%         cuts = NaN;
%     else
%         cuts = climit;
%         clear climit;
%         climit(1) = nanmin(cuts);
%         climit(2) = nanmax(cuts);
%     end
% end

%% prune data
if isreal(pixarr) ~= 1
    warning('Found complex data in pixarr. Taking real part.');
    pixarr = real(pixarr);
end

% fprintf(1,'In %s extrema are %g %g +/- %g \n,',mfilename,nanmin(nanmin(pixarr)),nanmax(nanmax(pixarr)),std(colvec(pixarr)));

%axis off;
switch idatatype
    case 0 %% values are wrapped phase in cycles on [-0.5 to +0.5]
        % scale to indices into a color table with 64 labels
        % add 33 to center in middle of color bar.
        % image (sans "s") does not use automatic scaling
        %image([xutmmin xutmmax],[yutmmax yutmmin],floor(33+64*pixarr)); %
        image([xutmmin xutmmax]/1000.,[yutmmax yutmmin]/1000.,floor(33+64*pixarr)); %
    case -1 %% values are gradients
        % use automatic scaling within limits
        imagesc([xutmmin xutmmax]/1000.,[yutmmax yutmmin]/1000.,pixarr,climits); % 
    case 2 %% values are (unwrapped) range change in meters
       % use automatic scaling within limits
        imagesc([xutmmin xutmmax]/1000.,[yutmmax yutmmin]/1000.,pixarr,climits);
    case 3 %% plot individual pixels as dots
        % assume that pixarr is actually an array with N rows and 3 columns
        E=colvec(pixarr(:,1)/1.e3);  % UTM easting in km
        N=colvec(pixarr(:,2)/1.e3);  % UTM northing in km
        V=colvec(pixarr(:,3));       % field values in (arbitrary) input units
        % try to guess pixel dimension in meters
        % dE=nanmedian(abs(diff(E)));
        % dN=nanmedian(abs(diff(N)));
        % find a number in points
        circleSize=marksize;
        circleSize=nanmin([ceil((xutmmax-xutmmin)/sqrt(numel(V))),ceil((yutmmax-yutmmin)/sqrt(numel(V)))]);
        circleSize=nanmin([circleSize,3]);

        if numel(climits) > 2
            % stretch color table
            [H, Vstretch] = contourfnu2(E,N,V,'levels',climits,'cmap',ctab,'method','scatter','symsize',circleSize,'symbol','s','pos_colorbar','none');
        else
            % use automatic scaling within limits
            scatter(E,N,circleSize,V,'filled'); % plot in km and mm
            caxis(climits);
        end
    otherwise
        error(sprintf('unknown idatatype %d\n',idatatype));
end

hold on;
axis xy; 
% 20200205 try to fix aspect ratio
axis equal;
axis tight;
% 2021/06/22  % draw a frame around
box on
set(gca,'Clipping','on','LineWidth',1);

if isa(ctab,'numeric')
    colormap(ctab)
end
if (drawxlabels == 0)
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel','');
    set(gca,'XColor',labelcolor);  % 2011-JUL-04
    %axis off;
else
    axis on;
    set(gca,'FontName','Helvetica','FontWeight','bold');
    xticklabels = get(gca,'XTickLabel');
    [nr,nc]=size(xticklabels);
    %xticklabels
    %whos xticklabels
    xticks = get(gca,'XTick');
    %whos xticks
    if nr*nc > 5
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',xticklabels(1:2:nr,:));
        %get(gca,'XTickLabel')
        set(gca,'XTick',xticks(1:2:nr));
        %get(gca,'XTick')
    end
    set(gca,'XColor',labelcolor);  % 2011-JUL-04
end

if (drawylabels == 0)
    set(gca,'YTickLabelMode','manual');
    set(gca,'YTickLabel','');
    set(gca,'YColor',labelcolor);  % 2011-JUL-04
    %axis off;
else
    %axis on;
    set(gca,'FontName','Helvetica','FontWeight','bold');
    yticklabels=get(gca,'YTickLabel');
    [nr,nc]=size(yticklabels);
    %yticklabels
    %whos yticklabels
    yticks = get(gca,'ytick');
    %whos yticks
    if nr*nc > 5
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',yticklabels(1:2:nr,:));
        %get(gca,'YTickLabel')
        set(gca,'YTick',yticks(1:2:nr));
        %get(gca,'YTick')
    end
    set(gca,'YColor',labelcolor);  % 2011-JUL-04
end

ha=text(0.02,0.98,cornerlabel...
,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
,'HorizontalAlignment','Left','VerticalAlignment','Top'...
,'BackgroundColor',[1 1 1 ]); % white box underneath

%ha=text(0.50,0.99,titlestr...
%,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
%,'HorizontalAlignment','Center','VerticalAlignment','Top'...
%,'BackgroundColor',[1 1 1 ]); % white box underneath

% ha=text(0.95,0.98,datelabel...
% ,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
% ,'HorizontalAlignment','Right','VerticalAlignment','Top'...
% ,'BackgroundColor',[1 1 1 ]); % white box underneath



% %% plot markers
% if marksize > 0 && numel(dotx) > 0 && numel(doty) > 0 && numel(mysym) == 3
%   plot(dotx/1000,doty/1000,mysym(1:2),'MarkerSize',marksize,'MarkerFaceColor',mysym(3)); 
% end

%% 20180529 add this line to avoid issue on unix systems with missing panel
%drawnow;

% % draw colorbar
% if cbar == 1
%     Hcb = colorbar;
%     Hcb.Label.String=datalabel;
% end

%% return with status
%istat = 0;
% return with handle to axis
H=gca;

% return graphics handle
if nargout >= 1
    varargout{1} = H;
end

% return values rescaled to stretch histogram
if nargout >= 2 && numel(isfinite(Vstretch)) > 0
    varargout{2} = Vstretch;
end


return
end
 

