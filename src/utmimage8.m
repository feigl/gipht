function h1 = utmimage8(im1,im2,im3,im4,im5,im6,im7,im8...
    ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
    ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes...
    ,idatatype,datalabel)
% function h1 = utmimage8(im1,im2,im3,im4,im5,im6,im7,im8...
%                        ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
%             ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes...
%             ,idatatype,datalabel)
%  1=a   5=e  obs on [-0.5,+0.5]
%  2=b   6=f  mod on [-0.5,+0.5]
%  3=c   7=g  res on [-0.5,+0.5]
%  4=d   8=h  costs on [0.0, 0.5]
% Last Modified 20160817

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
%rat = (8/11)*(yutmdif/xutmdif);
% 20180529 try to fix missing panel f
rat = (8.5/11.)*(yutmdif/xutmdif);
if isfinite(rat)
    if rat > 1.0
        warning(sprintf('Scale ratio (%10.4f) is greater than 1.0\n',rat));
        rs = min([1.0, 0.90/rat]); % rescale factor 2011-JUL-04
    elseif rat < 0.8
        warning(sprintf('Scale ratio (%10.4f) is less than 0.8\n',rat));
        %rs = min([1.0 1.5*rat]); % rescale factor 20161020
        rs = 1.0; % rescale factor 20180529
    else
        rs = 1.0;
    end
else
    rat = 1.0;
    rs  = 1.0;
end

nodot = 100;

h1=figure;
%% These are inches for US letter size paper
%set(h1,'PaperPosition', [0.25 0.25 8 8]);
%set(h1,'PaperPosition', [0.25 0.25 8 11]);
set(h1,'PaperPosition', [0 0 8.5 11]);

set(h1,'PaperPositionMode','manual');
%set(h1,'PaperSize',[8.5 8.5]);
set(h1,'PaperSize',[8.5 11.0]);

%set(h1,'Position',[1 1 1000 1000*rat]);
% 2021/06/22
set(h1,'Position',[1 1 1000 1000*11.0/8.5]);



%% parse levels and extrema
switch idatatype
    case 0
        labu = 'cycles';
        if abs(climit(2)-0.5) < 1e-6
            labt = '+1/2 (cycle)';
        else
            labt = sprintf('%+4.1f',climit(1));
        end
        if abs(climit(1)+0.5) < 1e-6
            labb = '-1/2 (cycle)';
        else
            labb = sprintf('%+4.1f',climit(2));
        end
    case -1
        labu = 'strain';
        zmin = nanmin(climit);
        zmax = nanmax(climit);
        labt = sprintf('%+5.1E',zmax);
        labb = sprintf('%+5.1E',zmin);
    case {2,3}
        %labu = 'mm';
        labu = datalabel;
        zmin = nanmin(climit);
        zmax = nanmax(climit);
        if abs(zmin) < 1000 && abs(zmin) < 1000 && abs(zmin) > 1 && abs(zmin) > 1
            labt = sprintf('%+5.0f',zmax);
            labb = sprintf('%+5.0f',zmin);
         elseif abs(zmin) < 1 && abs(zmin) < 1 && abs(zmin) > 0.1 && abs(zmin) > 0.1
            labt = sprintf('%+5.1f',zmax);
            labb = sprintf('%+5.1f',zmin);
        else
            labt = sprintf('%+5.1E',zmax);
            labb = sprintf('%+5.1E',zmin);
        end      
    otherwise
        error(sprintf('unknown idatatype %d\n',idatatype));
end


% fprintf(1,'In %s extrema are %g %g +/- %g \n,',mfilename,nanmin(nanmin(im1)),nanmax(nanmax(im1)),std(colvec(im1)));

colormap(ctab);

% draw paneles in smart order to see tick labels
% function istat = utmimage(pixarr,xutmmin,xutmmax,yutmmin,yutmmax...
%    ,titlestr,cornerlabel,climit,dotx,doty,ctab,mysym,marksize...
%    ,drawxlabels,drawylabels,datelabel...
%    ,idatatype,datalabel)

subplot('position',[0.325  0.775*rat  0.225 0.225*rat]*rs);%drawnow;
utmimage(im5,xutmmin,xutmmax,yutmmin,yutmmax,tl5,'e',climit,dotxutm, dotyutm,ctab,mysyms{5},marksizes(5),0,0,'',idatatype,datalabel);
title(tl5,'FontName','Helvetica','FontWeight','bold');
if cbar ==1
    ha=colorbar('Position',[0.550 0.775*rat 0.015 0.2250*rat]*rs,'YTickLabel',[]);
    ha=text(1.1,1.1, datalabel,'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor); 
    ha=text(1.1,0.05,   labb,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.95,   labt,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.50,  'Obs',  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor...
        ,'rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','Color',labelcolor);  
end

%% TODO - figure out why this panel disappears
subplot('position',[0.325  0.550*rat  0.225 0.225*rat]*rs);%drawnow;
utmimage(im6,xutmmin,xutmmax,yutmmin,yutmmax,tl6,'f',climit,dotxutm, dotyutm,ctab,mysyms{6},marksizes(6),0,0,'',idatatype,datalabel);
if cbar ==1
    ha=colorbar('Position',[0.550 0.550*rat 0.015 0.2250*rat]*rs,'YTickLabel',[]);
    ha=text(1.1,0.05,   labb,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.95,   labt,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.50,  'Mod',  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor...
        ,'rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','Color',labelcolor);  
end

subplot('position',[0.325  0.325*rat  0.225 0.225*rat]*rs);%drawnow;
utmimage(im7,xutmmin,xutmmax,yutmmin,yutmmax,tl7,'g',climit,dotxutm, dotyutm,ctab,mysyms{7},marksizes(7),0,0,'',idatatype,datalabel);
if cbar ==1
    ha=colorbar('Position',[0.550 0.325*rat 0.015 0.2250*rat]*rs,'YTickLabel',[]);
    ha=text(1.1,0.05,   labb,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.95,   labt,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.50,  'Res',  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor...
        ,'rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','Color',labelcolor);  
end

subplot('position',[0.325  0.1*rat  0.225 0.225*rat]*rs);%drawnow;
utmimage(im8,xutmmin,xutmmax,yutmmin,yutmmax,tl8,'h',climit,dotxutm, dotyutm,ctab,mysyms{8},marksizes(8),0,0,'',idatatype,datalabel);
ha=text(+0.2,-0.03 , 'Easting (km)', 'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','rotation', 0,'HorizontalAlignment','Left','VerticalAlignment','Top');
set(ha,'Color',labelcolor);  % 2011-JUL-04
if cbar == 1
    ha=colorbar('Position',[0.550 0.1*rat 0.015 0.2250*rat]*rs,'YTickLabel',[]);
    ha=text(1.1,0.50,   labb,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.95,   labt,  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor);  
    ha=text(1.1,0.50,  'Dev',  'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','Color',labelcolor...
        ,'rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','Color',labelcolor);  
end

subplot('position',[0.1  0.775*rat   0.225 0.225*rat]*rs);%drawnow;
utmimage(im1,xutmmin,xutmmax,yutmmin,yutmmax,tl1,'a',climit,dotxutm,dotyutm,ctab,mysyms{1},marksizes(1),0,0,'',idatatype,datalabel);
title(tl1,'FontName','Helvetica','FontWeight','bold');

subplot('position',[0.1  0.550*rat   0.225 0.225*rat]*rs);%drawnow;
utmimage(im2,xutmmin,xutmmax,yutmmin,yutmmax,tl2,'b',climit,dotxutm,dotyutm,ctab,mysyms{2},marksizes(2),0,0,'',idatatype,datalabel);

subplot('position',[0.1  0.325*rat   0.225 0.225*rat]*rs);%drawnow;
utmimage(im3,xutmmin,xutmmax,yutmmin,yutmmax,tl3,'c',climit,dotxutm,dotyutm,ctab,mysyms{3},marksizes(3),0,0,'',idatatype,datalabel);
ha=text(-0.05,+0.225, 'Northing (km)','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
set(ha,'Color',labelcolor);  % 2011-JUL-04

subplot('position',[0.1  0.1*rat   0.225 0.225*rat]*rs);%drawnow;
utmimage(im4,xutmmin,xutmmax,yutmmin,yutmmax,tl4,'d',climit,dotxutm, dotyutm,ctab,mysyms{4},marksizes(4),1,1,'',idatatype,datalabel);

subplot('position',[0.325  0.1*rat   0.225 0.225*rat]*rs);%drawnow;
utmimage(im8,xutmmin,xutmmax,yutmmin,yutmmax,tl8,'h',climit,dotxutm, dotyutm,ctab,mysyms{8},marksizes(8),1,1,'',idatatype,datalabel);

%% label with title at top
%subplot('position',[0.1 1.0*rat 0.5 0.05*rat]*rs,'Units','normalized');%drawnow;
%% absolute top
subplot('position',[0.1 0.95 0.5 0.05],'Units','normalized');%drawnow;
axis off
% coordinates for text are inside the rectangle defined by subplot above
text(0.1,0.1,titlestr ...
    ,'FontName','Helvetica','Fontsize',10,'FontWeight','Bold' ...
    ,'HorizontalAlignment','Left','VerticalAlignment','Bottom' ...
    ,'Clipping','off' ...
    ,'Units','Normalized','rotation', 0 ...
    ,'margin',2 ...
    ,'Interpreter','none');
drawnow;

return;
end





