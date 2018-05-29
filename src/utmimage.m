function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax...
   ,titlestr,cornerlabel,climit,dotx,doty,ctab,mysym,marksize...
   ,drawxlabels,drawylabels,datelabel...
   ,idatatype,datalabel)
% function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax...
%    ,titlestr,cornerlabel,climit,dotx,doty,ctab,mysym,marksize...
%    ,drawxlabels,drawylabels,datelabel...
%    ,idatatype,datalabel)


% initialize to return an error
istat = -1;

% if nargin < 8 || exist('climit','var') == 0 || exist('dotx','var') == 0 || exist('doty','var') == 0
%    climit  = [-0.5 0.5];
%    dotx = 0;
%    doty = 0;
% end
% if nargin < 11 || exist('ctab','var') == 0
%    ctab = 'jet';
% end
% if nargin < 12 || exist('mysym','var') == 0
%    mysym = 'ko';  % default is draw a black circle, but do not connect
% end
% if nargin < 13 || exist('marksize','var') == 0
%    marksize = 5;
% end
% if nargin < 14 || exist('drawxlabels','var') == 0
%    drawxlabels=0;
% end
% if nargin < 15 || exist('drawylabels','var') == 0
%    drawylabels=0;
% end
% if nargin < 16 || exist('datelabel','var') == 0
%    datelabel='';
% end
% if nargin < 17 || exist('idatatype','var') == 0
%    idatatype = 0; % wrapped phase is default
% end

%labelcolor = [0.5 0.5 0.5]; % gray
labelcolor = [0.0 0.0 0.0]; % black

% UTM coordinates in km
xutmmin = xutmmin/1000;
xutmmax = xutmmax/1000;
yutmmin = yutmmin/1000;
yutmmax = yutmmax/1000;

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


% prune
if isreal(pixarr) ~= 1
    warning('Found complex data in pixarr. Taking real part.');
    pixarr = real(pixarr);
end

% fprintf(1,'In %s extrema are %g %g +/- %g \n,',mfilename,nanmin(nanmin(pixarr)),nanmax(nanmax(pixarr)),std(colvec(pixarr)));

switch idatatype
    case 2 %% values are (unwrapped) range change in meters          
       %imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr); % use automatic scaling
       imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr,climit); % use automatic scaling within limits
    case 0 %% values are wrapped phase in cycles on [-0.5 to +0.5]
           % scale to indices into a color table with 64 labels 
           % add 33 to center in middle of color bar. 
           % image (sans "s") does not use automatic scaling
       image([xutmmin xutmmax],[yutmmax yutmmin],floor(33+64*pixarr)); %
    case -1 %% values are gradients
       %imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr); % use automatic scaling
       imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr,climit); % use automatic scaling within limits
    otherwise
        error(sprintf('unknown idatatype %d\n',idatatype));
 end

if (drawxlabels == 0)
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel','');
else
    set(gca,'FontName','Helvetica','FontWeight','bold');
    xticklabels = get(gca,'XTickLabel');
    [nr,nc]=size(xticklabels);
    %xticklabels
    %whos xticklabels
    xticks = get(gca,'XTick');
    %whos xticks
    if nr*nc > 10
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',xticklabels(1:2:nr,:));
        %get(gca,'XTickLabel')
        set(gca,'XTick',xticks(1:2:nr));
        %get(gca,'XTick')
    end    
end
set(gca,'XColor',labelcolor);  % 2011-JUL-04

if (drawylabels == 0)
    set(gca,'YTickLabelMode','manual');
    set(gca,'YTickLabel','');
else
    set(gca,'FontName','Helvetica','FontWeight','bold');
    yticklabels=get(gca,'YTickLabel');
    [nr,nc]=size(yticklabels);
    %yticklabels
    %whos yticklabels
    yticks = get(gca,'ytick');
    %whos yticks
    if nr*nc > 10
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',yticklabels(1:2:nr,:));
        %get(gca,'YTickLabel')
        set(gca,'YTick',yticks(1:2:nr));
        %get(gca,'YTick')
    end
    
end
set(gca,'YColor',labelcolor);  % 2011-JUL-04

ha=text(0.02,0.98,cornerlabel...
,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
,'HorizontalAlignment','Left','VerticalAlignment','Top'...
,'BackgroundColor',[1 1 1 ]); % white box underneath

%ha=text(0.50,0.99,titlestr...
%,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
%,'HorizontalAlignment','Center','VerticalAlignment','Top'...
%,'BackgroundColor',[1 1 1 ]); % white box underneath

ha=text(0.95,0.98,datelabel...
,'Units','normalized','Clipping','off','FontName','Helvetica','margin',1,'FontWeight','bold'...
,'HorizontalAlignment','Right','VerticalAlignment','Top'...
,'BackgroundColor',[1 1 1 ]); % white box underneath

hold on;
axis xy; 

%% plot markers
if marksize > 0 && numel(dotx) > 0 && numel(doty) > 0
  plot(dotx/1000,doty/1000,mysym,'MarkerSize',marksize); hold on;
end

%% 20180529 add this line to avoid issue on unix systems with missing panel
drawnow;

%% return with status
istat = 0;


return
 

