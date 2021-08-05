function H1 = utmimage8landscape(im1,im2,im3,im4,im5,im6,im7,im8...
            ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
            ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes...
            ,idatatype, datalabel)
% function h1 = utmimage8landscape(im1,im2,im3,im4,im5,im6,im7,im8...
%                        ,tl1,tl2,tl3,tl4,tl5,tl6,tl7,tl8....
%             ,wesn,titlestr,climit,dotxutm,dotyutm,ctab,cbar,mysyms,marksizes...
%             ,'',idatatype, datalabel)
%
%          obs   mod  res  dev
%  initial 1=a   2=b  3=c  4=d
%  final   5=e   6=f  7=g  8=h
%
% Last Modified 2021/06/22 
% Tried unsuccessfully to remove white space between plots in horizontal
% dimension


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
%rat = (8.5/11)*(yutmdif/xutmdif);
% 2021/06/22
rat = (yutmdif/xutmdif)

if isfinite(rat)
    if rat > 1.0
        warning(sprintf('Scale ratio (%10.4f) is greater than 1.0\n',rat));
        rs = min([1.0, 0.90/rat]); % rescale factor 2011-JUL-04
    elseif rat < 0.8
        warning(sprintf('Scale ratio (%10.4f) is less than 0.8\n',rat));
        %rs = min([1.0 1.5*rat]); % rescale factor 20161020
        rs = 1.0; % rescale factor 20180529
    else
        rs = 0.8;  % fit on page
    end
else
    rat = 1.0;
    rs  = 1.0;
end

nodot = 100;

H1=figure;
orient(H1,'landscape');

%% make landscape -- no effect if printpdf is called afterwards
set(H1,'PaperUnits','inch');
set(H1,'PaperPositionMode','manual');
set(H1,'PaperPosition', [0. 0. 11 8.5]);
set(H1,'PaperOrientation', 'landscape');
%set(h1,'PaperSize',[8.5 11.0]);
set(H1,'PaperSize',[11.0 8.5]);
%set(h1,'Position',[1 1 1000 1000*rat]);
% 2021/05/22
set(H1,'Position',[1 1 1000*11./8.5 1000]);


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
        labt = sprintf('%+5.1E',zmin);
        labb = sprintf('%+5.1E',zmax);
    case {2,3}
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
labu = datalabel;

colormap(ctab);

% draw panels in smart order to see tick labels

subplot('Position',[0.825  0.325*rat   0.225 0.225*rat]*rs); 
%subplot('Position',[0.800  0.325*rat   0.225 0.225*rat]*rs); axis off; drawnow;
H=utmimage(im4,xutmmin,xutmmax,yutmmin,yutmmax,tl4,'d',climit,dotxutm, dotyutm,ctab,mysyms{4},marksizes(4),0,0,'',idatatype,datalabel);
% subplot(H,'Position',[0.800  0.325*rat   0.225 0.225*rat]*rs); axis off; drawnow;

if cbar == 1  
   ha=text(1.25 ,-0.850,labb            ,'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Right'  ,'VerticalAlignment','Bottom','rotation', 0);
   ha=text(1.25 , 0.000,labu            ,'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',90);    
   ha=text(1.25 , 0.850,labt            ,'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Right ' ,'VerticalAlignment','Top','rotation', 0); 
   ha=text(0.50 , 1.0      ,'Deviation' ,'Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.600  0.325*rat   0.225 0.225*rat]*rs);
%subplot('position',[0.625  0.325*rat   0.225 0.225*rat]*rs);
%pos=get(gca,'InnerPosition')
% subplot('position',[0.825-1.01*pos(3)  0.325*rat   0.225 0.225*rat]*rs); %drawnow;
%subplot('position',[0.600-0.015  0.325*rat   0.225 0.225*rat]*rs,'align'); axis off; drawnow;

% %https://www.mathworks.com/matlabcentral/answers/93012-how-do-i-decrease-the-margins-around-the-subplots-in-my-figure-in-matlab
% ax = gca;
% pos = ax.Position
% innerpos = ax.InnerPosition
% outerpos = ax.OuterPosition
% ti=outerpos-innerpos
%ti = ax.TightInset 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
%subplot('position',[0.600+ti(3), 0.325*rat,  0.225, 0.225*rat]*rs); axis off; drawnow;

H=utmimage(im3,xutmmin,xutmmax,yutmmin,yutmmax,tl3,'c',climit, dotxutm,dotyutm,ctab,mysyms{3},marksizes(3),0,0,'',idatatype,datalabel);
if cbar == 1  
   ha=text(0.5  , 1.0      ,'Residual','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.375  0.325*rat   0.225 0.225*rat]*rs); 
utmimage(im2,xutmmin,xutmmax,yutmmin,yutmmax,tl2,'b',climit, dotxutm,dotyutm,ctab,mysyms{2},marksizes(2),0,0,'',idatatype,datalabel);
if cbar == 1  
   ha=text(0.5  , 1.0      ,'Modeled','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
end

subplot('position',[0.150    0.325*rat   0.225 0.225*rat]*rs); 
H=utmimage(im1,xutmmin,xutmmax,yutmmin,yutmmax,tl1,'a',climit, dotxutm,dotyutm,ctab,mysyms{1},marksizes(1),0,0,'',idatatype,datalabel);
if cbar == 1  
   ha=text(0.50  ,  1.0      ,'Observed','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','Center' ,'VerticalAlignment','Bottom','rotation',0);    
   ha=text(-0.15,+0.1, 'Northing (km)','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Top');
end

subplot('position',[0.825  0.1*rat  0.225 0.225*rat]*rs); 
H=utmimage(im8,xutmmin,xutmmax,yutmmin,yutmmax,tl8,'h',climit,dotxutm, dotyutm,ctab,mysyms{8},marksizes(8),0,0,'',idatatype,datalabel);

subplot('position',[0.600  0.1*rat  0.225 0.225*rat]*rs);
H=utmimage(im7,xutmmin,xutmmax,yutmmin,yutmmax,tl7,'g',climit,dotxutm, dotyutm,ctab,mysyms{7},marksizes(7),0,0,'',idatatype,datalabel);

subplot('position',[0.375  0.1*rat  0.225 0.225*rat]*rs); 
H=utmimage(im6,xutmmin,xutmmax,yutmmin,yutmmax,tl6,'f',climit,dotxutm, dotyutm,ctab,mysyms{6},marksizes(6),0,0,'',idatatype,datalabel);

if cbar == 1
   ha=text(+0.2,-0.03    , 'Easting (km)','Units','normalized','Clipping','off','FontName','Helvetica','FontWeight','Bold','rotation', 0,'HorizontalAlignment','Left','VerticalAlignment','Top');
end

subplot('position',[0.150  0.1*rat  0.225 0.225*rat]*rs); 
H=utmimage(im5,xutmmin,xutmmax,yutmmin,yutmmax,tl5,'e',climit,  dotxutm, dotyutm,ctab,mysyms{5},marksizes(5),1,1,'',idatatype,datalabel);

% 2021/06/28 draw color bar last
if cbar == 1
   ha=colorbar('Position',[1.100 0.1*rat 0.015 0.450*rat]*rs); 
   set(ha,'YTickLabel',[]);
   drawnow;
end

% Wasted a lot of time trying unsuccessfully to make all this stuff work
%set(ha,'YLimMode','manual');
%set(ha,'YLim',[0 1]);
%set(ha,'YTickLabelMode','manual');
%set(ha,'YTickMode','manual');
%set(ha,'YTick',[0 0.1 0.2]);
%set(ha,'ZTickMode','manual','ZTick',[ybot+0.01],'ZTickLabelMode','manual','ZTickLabel','0.0');


% % leave this spot free for yutmer labelling
% % subplot('position',[0.9  0.1  0.1 0.9]);
% % axis off
% 
% subplot('position',[0.1 0.6*rat 0.1 0.02])
% axis off
% text(0.0,0.0,strrep(titlestr,'\n',' '),...
%    'FontName','Helvetica','Fontsize',10,'FontWeight','normal',...
%    'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
%    'Units','Normalized','rotation', 0,'BackgroundColor',[1 1 1],...
%    'margin',2);

%% label with title at top
subplot('position',[0.1 1.0*rat 0.5 0.05*rat]*rs,'Units','normalized'); axis off; drawnow;
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

% return current figure handle
H1=gcf;
return;
end





