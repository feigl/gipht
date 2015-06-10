function h1 = utmimage8landscape(im1,im2,im3,im4,im5,im6,im7,im8...
                       ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
            ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes)
% function h1 = utmimage8landscape(im1,im2,im3,im4,im5,im6,im7,im8...
%                        ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
%             ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes)
%
%          obs   mod  res  dev
%  initial 1=a   2=b  3=c  4=d
%  final   5=e   6=f  7=g  8=h
%
% Last Modified 2010-JUL-21


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
rat = (11/8.5)*(yutmdif/xutmdif);

if rat > 1.0
    warning(sprintf('Scale ratio (%10.4f) is greater than 1.0\n',rat));
    rs = min([1.0, 0.90/rat]); % rescale factor 2011-JUL-04
else
    rs = 1.0;
end

nodot = 100;

h1=figure;
%set(h1,'PaperPosition', [0.25 0.25 8 8]);
set(h1,'PaperPosition', [0. 0. 11 8.5]);

set(h1,'PaperPositionMode','manual');
%set(h1,'PaperSize',[8.5 8.5]);
set(h1,'PaperSize',[11.0 8.5]);


set(h1,'Position',[1 1 1000 1000*rat]);
%set(h1,'Position',[1 1 1000 1000]);

% parse levels and extrema
nlevels = numel(climit);
if nlevels == 2
    if abs(climit(2)-0.5) < 1e-6
        labt = '+1/2';
        labm = 'phase (cycle)';
    else
        labt = sprintf('%+4.1f',climit(2));
        labm = 'gradient';
    end
    if abs(climit(1)+0.5) < 1e-6
        labb = '-1/2';
    else
        labb = sprintf('%+4.1f',climit(1));
    end
else    
    labb = sprintf('%+5.0G',min(climit));
    labt = sprintf('%+5.0G',max(climit));
    labm = 'gradient';
end


colormap(ctab);


% draw panels in smart order to see tick labels

subplot('position',[0.825  0.325*rat   0.225 0.225*rat]*rs);
utmimage(im4,xutmmin,xutmmax,yutmmin,yutmmax,tl4,'d',climit,dotxutm, dotyutm,ctab,mysyms{4},marksizes(4),0,0);
if cbar == 1  
   ha=text(1.45 ,-0.850,labb            ,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Right'  ,'VerticalAlignment','Bottom','rotation', 0);
   ha=text(1.45 , 0.000,labm            ,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',90);    
   ha=text(1.45 , 0.850,labt            ,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Right ' ,'VerticalAlignment','Top','rotation', 0); 
   ha=text(0.50 , 1.0      ,'Deviation' ,'Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.600  0.325*rat   0.225 0.225*rat]*rs);
utmimage(im3,xutmmin,xutmmax,yutmmin,yutmmax,tl3,'c',climit, dotxutm,dotyutm,ctab,mysyms{3},marksizes(3),0,0);
if cbar == 1  
   ha=text(0.5  , 1.0      ,'Residual','Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.375  0.325*rat   0.225 0.225*rat]*rs);
utmimage(im2,xutmmin,xutmmax,yutmmin,yutmmax,tl2,'b',climit, dotxutm,dotyutm,ctab,mysyms{2},marksizes(2),0,0);
if cbar == 1  
   ha=text(0.5  , 1.0      ,'Modeled','Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.150    0.325*rat   0.225 0.225*rat]*rs);
utmimage(im1,xutmmin,xutmmax,yutmmin,yutmmax,tl1,'a',climit, dotxutm,dotyutm,ctab,mysyms{1},marksizes(1),0,0);
if cbar == 1  
   ha=text(0.50  ,  1.0      ,'Observed','Units','normalized','Clipping','off','FontName','Helvetica-Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
   ha=text(-0.25,+0.1, 'Northing (km)','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Top');
end

if cbar == 1
   ha=colorbar('Position',[1.050 0.1*rat 0.015 0.450*rat]*rs);
   set(ha,'YTickLabel',[]);
end

subplot('position',[0.825  0.1*rat  0.225 0.225*rat]*rs);%colormap('jet');
utmimage(im8,xutmmin,xutmmax,yutmmin,yutmmax,tl8,'h',climit,dotxutm, dotyutm,ctab,mysyms{8},marksizes(8),0,0);


subplot('position',[0.600  0.1*rat  0.225 0.225*rat]*rs);
utmimage(im7,xutmmin,xutmmax,yutmmin,yutmmax,tl7,'g',climit ,dotxutm, dotyutm,ctab,mysyms{7},marksizes(7),0,0);


subplot('position',[0.375  0.1*rat  0.225 0.225*rat]*rs);
utmimage(im6,xutmmin,xutmmax,yutmmin,yutmmax,tl6,'f',climit  ,dotxutm, dotyutm,ctab,mysyms{6},marksizes(6),0,0);

if cbar == 1
   ha=text(+0.2,-0.03    , 'Easting (km)','Units','normalized','Clipping','off','FontName','Helvetica-Bold','rotation', 0,'HorizontalAlignment','Left','VerticalAlignment','Top');
end

subplot('position',[0.15  0.1*rat  0.225 0.225*rat]*rs);
utmimage(im5,xutmmin,xutmmax,yutmmin,yutmmax,tl5,'e',climit,  dotxutm, dotyutm,ctab,mysyms{5},marksizes(5),1,1);


% Wasted a lot of time trying unsuccessfully to make all this stuff work
%set(ha,'YLimMode','manual');
%set(ha,'YLim',[0 1]);
%set(ha,'YTickLabelMode','manual');
%set(ha,'YTickMode','manual');
%set(ha,'YTick',[0 0.1 0.2]);
%set(ha,'ZTickMode','manual','ZTick',[ybot+0.01],'ZTickLabelMode','manual','ZTickLabel','0.0');


% leave this spot free for yutmer labelling
% subplot('position',[0.9  0.1  0.1 0.9]);
% axis off

subplot('position',[0.1 0.6*rat 0.1 0.02])
axis off
text(0.0,0.0,strrep(titlestr,'\n',' '),...
   'FontName','Helvetica','Fontsize',10,'FontWeight','normal',...
   'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
   'Units','Normalized','rotation', 0,'BackgroundColor',[1 1 1],...
   'margin',2);

return;





