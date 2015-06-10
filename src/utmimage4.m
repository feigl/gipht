function h1 = utmimage4(im1,im2,im3,im4 ...
                       ,tl1,tl2,tl3,tl4 ....
            ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes)
%  1=a     obs on [-0.5,+0.5]
%  2=b     mod on [-0.5,+0.5]
%  3=c     res on [-0.5,+0.5]
%  4=d     costs on [0.0, 0.5]
% Last Modified 2011-JUL-13

%labelcolor = [0.5 0.5 0.5]; % Gray
labelcolor = [0. 0. 0.]; % Black

% if nargin < 19
%    cbar = 0;
% end
% if nargin < 20
%    for i=1:8
%       mysyms{i} = 'k-';
%    end
% end
% if nargin < 21
%    for i=1:8
%       marksizes(i) = 1;
%    end
% end
xutmmin = wesn(1);
xutmmax = wesn(2);
yutmmin = wesn(3);
yutmmax = wesn(4);

yutmmid = (yutmmax+yutmmin)/2.0;
yutmdif = yutmmax-yutmmin;
xutmdif = xutmmax-xutmmin;
rat = (8/11)*(yutmdif/xutmdif);
if rat > 1.0
    warning(sprintf('Scale ratio (%10.4f) is greater than 1.0\n',rat));
    rs = min([1.0, 0.90/rat]); % rescale factor 2011-JUL-04 
else
    rs = 1.0;
end
    
nodot = 100;

h1=figure;
%set(h1,'PaperPosition', [0.25 0.25 8 8]);
set(h1,'PaperPosition', [0.25 0.25 8 11]);

set(h1,'PaperPositionMode','manual');
%set(h1,'PaperSize',[8.5 8.5]);
set(h1,'PaperSize',[8.5 11.0]);

set(h1,'Position',[1 1 1000 1000*rat]);


% parse levels and extrema
if abs(nanmax(climit)-nanmin(climit)) > 1.0
    labu = 'mm';
    labt = sprintf('%4.0f',max(climit));
    labb = sprintf('%4.0f',min(climit));
else
    labu = 'cycles';
    labt = '+1/2';
    labb = '-1/2';   
end
colormap(ctab);

% draw panels in smart order to see tick labels

subplot('position',[0.325  0.775*rat  0.225 0.225*rat]*rs);colormap(ctab);
utmimage(im2,xutmmin,xutmmax,yutmmin,yutmmax,tl2,'b',climit,dotxutm, ...
         dotyutm,ctab,mysyms{2},marksizes(2),0,0);
title(tl2,'FontName','Helvetica-Bold');
if cbar ==1
    ha=colorbar('Position',[0.550 0.775*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.05,labb,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.50,labu,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end


subplot('position',[0.325  0.550*rat  0.225 0.225*rat]*rs);
utmimage(im4,xutmmin,xutmmax,yutmmin,yutmmax,tl4,'d',climit,dotxutm, dotyutm,ctab,mysyms{4},marksizes(4),0,0);
if cbar ==1
    ha=colorbar('Position',[0.550 0.550*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.05,labb,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.50, labu,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(+0.2,-0.05*rat , 'Easting (km)', 'Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation', 0,'HorizontalAlignment','Left','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04    
end

subplot('position',[0.1  0.775*rat   0.225 0.225*rat]*rs);
utmimage(im1,xutmmin,xutmmax,yutmmin,yutmmax,tl1,'a',climit,dotxutm,dotyutm,ctab,mysyms{1},marksizes(1),0,0);
title(tl1,'FontName','Helvetica-Bold');
if cbar == 1
    ha=text(-0.05,+0.375*rat, 'Northing (km)','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.1  0.550*rat   0.225 0.225*rat]*rs);
utmimage(im3,xutmmin,xutmmax,yutmmin,yutmmax,tl3,'c',climit,dotxutm,dotyutm,ctab,mysyms{3},marksizes(3),1,1);


subplot('position',[0.65 0.1 0.02 0.9])
axis off
text(0.0,0.0,titlestr,...
    'FontName','Helvetica','Fontsize',10,'FontWeight','normal',...
    'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
    'Units','Normalized','rotation', 90,'BackgroundColor',[1 1 1],...
    'margin',2);


return;





