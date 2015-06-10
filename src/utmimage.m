function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax,...
   titlestr,cornerlabel,climit,dotx,doty,ctab,mysym,marksize,...
                          drawxlabels,drawylabels,datelabel)
% function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax,...
%    titlestr,cornerlabel,climit,dotx,doty,ctab,mysym,marksize)
% make a map of an image
% Last Modified:
% 2009-JUN-18 - fix labels
% 2011-NOV-21 - swap order of ymin ymax for MATLAB release 2011b
% 2012-OCT-25 - add datelabel argument

% initialize to return an error
istat = -1;

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

%labelcolor = [0.5 0.5 0.5]; % gray
labelcolor = [0.0 0.0 0.0]; % black

% UTM coordinates in km
xutmmin = xutmmin/1000;
xutmmax = xutmmax/1000;
yutmmin = yutmmin/1000;
yutmmax = yutmmax/1000;

% make this routine robust
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

% parse levels and extrema
iok = find(isfinite(climit)==1);
climit = climit(iok);
nlevels = numel(climit);
if nlevels < 2
   cuts = NaN;
   climit(1) = -0.5;
   climit(2) = +0.5;   
elseif nlevels == 2
   cuts = NaN;
else 
   cuts = climit;
   clear climit;
   climit(1) = nanmin(cuts);
   climit(2) = nanmax(cuts);
end

% prune
if isreal(pixarr) ~= 1
    warning('Found complex data in pixarr. Taking real part.');
    pixarr = real(pixarr);
end

%KF 2009-MAR-29 imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr); 
%imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr,climit); 
pmin=floor(256*min(min(pixarr)));
pmax= ceil(256*max(max(pixarr)));
if pmin == 0 && pmax == 128
%     fprintf(1,'Plotting as cost values ranging from 0 to 127 without scaling.\n');
    image([xutmmin xutmmax],[yutmmax yutmmin],floor(33+64*pixarr)); % values are costs [0, 127]
    %image([xutmmin xutmmax],[yutmmin yutmmax],floor(33+64*pixarr)); axis xy
elseif pmin == -128 && pmax == 127
%     fprintf(1,'Plotting as phase values ranging from -128 to +127 without scaling.\n');    
    image([xutmmin xutmmax],[yutmmax yutmmin],floor(33+64*pixarr)); % values are phases [-127, +128]
   %image([xutmmin xutmmax],[yutmmin yutmmax],floor(33+64*pixarr)); axis xy % values are phases [-127, +128]
   
% elseif abs(climit(1)+0.5) > 1.0E-3 || abs(climit(2)-0.5) > 1.0E-3
% elseif nlevels > 2
%     fprintf(1,'%s equalizing histogram.\n',mfilename);
%     pixeqh = equalize_hist(pixarr,64);
% %     maxabs = max(max(abs(pixeqh)));
% %     climit = [-maxabs,maxabs];
%     imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixeqh,climit);
elseif nlevels > 2
%     fprintf(1,'%s indexing\n',mfilename);
   image([xutmmin xutmmax],[yutmmax yutmmin],pixarr);
   %image([xutmmin xutmmax],[yutmmin yutmmax],pixarr);axis xy;
else
    %      fprintf(1,'%s scaling\n',mfilename);
    %      imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr);
    %       fprintf(1,'%s scaling with climit %f %f\n',mfilename,min(climit),max(climit));
    try
        imagesc([xutmmin xutmmax],[yutmmax yutmmin],pixarr,climit);
        %imagesc([xutmmin xutmmax],[yutmmin yutmmax],pixarr,climit); axis xy;
    catch
        warning('imagesc failed');
        istat = -1;
        return
    end
end

% %equalize histogram
% if pselect == 7
%     size(imE)
%     [ctab, imE] = histeq(imE,64,1,'jet');
%     size(imE)
% end

if (drawxlabels == 0)
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel','');
else
    set(gca,'FontName','Helvetica-Bold');
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
    set(gca,'FontName','Helvetica-Bold');
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

ha=text(0.01,0.99,cornerlabel...
,'Units','normalized','Clipping','off','FontName','Helvetica-Bold'...
,'HorizontalAlignment','Left','VerticalAlignment','Top'...
,'BackgroundColor',[1 1 1 ]); % white box underneath

%ha=text(0.50,0.99,titlestr...
%,'Units','normalized','Clipping','off','FontName','Helvetica-Bold'...
%,'HorizontalAlignment','Center','VerticalAlignment','Top'...
%,'BackgroundColor',[1 1 1 ]); % white box underneath

ha=text(0.95,0.99,datelabel...
,'Units','normalized','Clipping','off','FontName','Helvetica-Bold'...
,'HorizontalAlignment','Right','VerticalAlignment','Top'...
,'BackgroundColor',[1 1 1 ]); % white box underneath


hold on;
axis xy; 

%colormap(ctab)

% if doty(numel(doty)) > 90
%    axis off;
% end

% plot markers, with coordinates in KILOMETERS
if marksize > 0 && numel(dotx) > 0 && numel(doty) > 0
  %plot(dotxutm,dotyutm,'ko','MarkerSize',5); hold on;
  plot(dotx,doty,mysym,'MarkerSize',marksize); hold on;
end



istat = 0;

return
 

