function h1 = utmimage8(im1,im2,im3,im4,im5,im6,im7,im8...
                       ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
            ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes)
% function h1 = utmimage8(im1,im2,im3,im4,im5,im6,im7...
%                        ,tl1,tl2,tl3,tl4,tl5,tl6,im8....
%             ,wesn,titlestr,climit,dotxutm,nodot,ctab,cbar)
%
%  1=a   5=e  obs on [-0.5,+0.5]
%  2=b   6=f  mod on [-0.5,+0.5]
%  3=c   7=g  res on [-0.5,+0.5]
%  4=d   8=h  costs on [0.0, 0.5]
% Last Modified 2011-JUL-04

%labelcolor = [0.5 0.5 0.5]; % Gray
labelcolor = [0. 0. 0.]; % Black

if nargin < 19
   cbar = 0;
end
if nargin < 20
   for i=1:8
      mysyms{i} = 'k-';
   end
end
if nargin < 21
   for i=1:8
      marksizes(i) = 1;
   end
end
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
nlevels = numel(climit);
% parse levels and extrema
nlevels = numel(climit);
if nlevels == 2
    if abs(climit(2)-0.5) < 1e-6
        labt = '+1/2 (cycle)';
    else
        labt = sprintf('%+4.1f',climit(2));
    end
    if abs(climit(1)+0.5) < 1e-6
        labb = '-1/2 (cycle)';
    else
        labb = sprintf('%+4.1f',climit(1));
    end
else    
    labb = sprintf('%+5.0G',min(climit));
    labt = sprintf('%+5.0G',max(climit));
end


colormap(ctab);

% draw paneles in smart order to see tick labels

subplot('position',[0.325  0.775*rat  0.225 0.225*rat]*rs);colormap(ctab);
utmimage(im5,xutmmin,xutmmax,yutmmin,yutmmax,tl5,'e',climit,dotxutm, ...
         dotyutm,ctab,mysyms{5},marksizes(5),0,0);
title(tl5,'FontName','Helvetica-Bold');
if cbar ==1
    ha=colorbar('Position',[0.550 0.775*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.05,labb,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.50, 'Obs','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.325  0.550*rat  0.225 0.225*rat]*rs);
utmimage(im6,xutmmin,xutmmax,yutmmin,yutmmax,tl6,'f',climit,dotxutm, dotyutm,ctab,mysyms{6},marksizes(6),0,0);
if cbar ==1
    ha=colorbar('Position',[0.550 0.550*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.05,labb,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.50, 'Mod','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.325  0.325*rat  0.225 0.225*rat]*rs);
utmimage(im7,xutmmin,xutmmax,yutmmin,yutmmax,tl7,'g',climit,dotxutm, dotyutm,ctab,mysyms{7},marksizes(7),0,0);
if cbar ==1
    ha=colorbar('Position',[0.550 0.325*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.05,labb,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.50, 'Res','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.325  0.1*rat  0.225 0.225*rat]*rs);%colormap('jet');
utmimage(im8,xutmmin,xutmmax,yutmmin,yutmmax,tl8,'h',climit,dotxutm, dotyutm,ctab,mysyms{8},marksizes(8),0,0);
if cbar == 1
    ha=colorbar('Position',[0.550 0.1*rat 0.015 0.2250*rat]*rs);
    set(ha,'YTickLabel',[]);
    ha=text(1.1,0.50,   '0','Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.05, 'Dev','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(1.1,0.95,labt,'Units','normalized','Clipping','off','FontName','Helvetica-Bold');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
    ha=text(+0.2,-0.03 , 'Easting (km)', 'Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation', 0,'HorizontalAlignment','Left','VerticalAlignment','Top');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.1  0.775*rat   0.225 0.225*rat]*rs);
utmimage(im1,xutmmin,xutmmax,yutmmin,yutmmax,tl1,'a',climit,dotxutm,dotyutm,ctab,mysyms{1},marksizes(1),0,0);
title(tl1,'FontName','Helvetica-Bold');

subplot('position',[0.1  0.550*rat   0.225 0.225*rat]*rs);
utmimage(im2,xutmmin,xutmmax,yutmmin,yutmmax,tl2,'b',climit,dotxutm,dotyutm,ctab,mysyms{2},marksizes(2),0,0);

subplot('position',[0.1  0.325*rat   0.225 0.225*rat]*rs);
utmimage(im3,xutmmin,xutmmax,yutmmin,yutmmax,tl3,'c',climit,dotxutm,dotyutm,ctab,mysyms{3},marksizes(3),0,0);
if cbar == 1
    ha=text(-0.1,+0.225, 'Northing (km)','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
    set(ha,'Color',labelcolor);  % 2011-JUL-04
end

subplot('position',[0.1  0.1*rat   0.225 0.225*rat]*rs);%colormap('jet');
utmimage(im4,xutmmin,xutmmax,yutmmin,yutmmax,tl4,'d',climit,dotxutm, dotyutm,ctab,mysyms{4},marksizes(4),1,1);



% Wasted a lot of time trying unsuccessfull to make all this stuff work
%set(ha,'YLimMode','manual');
%set(ha,'YLim',[0 1]);
%set(ha,'YTickLabelMode','manual');
%set(ha,'YTickMode','manual');
%set(ha,'YTick',[0 0.1 0.2]);
%set(ha,'ZTickMode','manual','ZTick',[ybot+0.01],'ZTickLabelMode','manual','ZTickLabel','0.0');




% leave this spot free for yutmer labelling
% subplot('position',[0.9  0.1  0.1 0.9]);
% axis off




subplot('position',[0.8 0.1 0.02 0.9])
axis off
text(0.0,0.0,strrep(titlestr,'\n',' '),...
    'FontName','Helvetica','Fontsize',10,'FontWeight','normal',...
    'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
    'Units','Normalized','rotation', 90,'BackgroundColor',[1 1 1],...
    'margin',2);


return;





