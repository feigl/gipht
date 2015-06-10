function imout=plot_patch_quads(qqq,qi1,qi2,qj1,qj2,cuts,cmap,tstr,zstr,drawcolorbar,nrows,ncols,drawit)
%function imout=plot_patch_quads(qqq,qi1,qi2,qj1,qj2,cuts,cmap,tstr,zstr,drawcolorbar,drawit)
% make a figure with patches colored based on cut color table
% return imout with indices into colormap
% draw figure if drawcolorbar == 1
% I indices refer to ROWS
% J indices refer to COLS
% Kurt Feigl 2010-JUL-08
% Revised 2014-JAN-07
% I indices refer to COLUMNS
% J indices refer to ROWS

if exist('drawit','var') == 0
    drawit = 0;
end

if drawit == 1
    h=figure;hold on;
end
iimin=nanmin(nanmin(qi1));
iimax=nanmax(nanmax(qi2));
jjmin=nanmin(nanmin(qj1));
jjmax=nanmax(nanmax(qj2));

% nr=ceil(abs(jmax-jmin)+1);
% nc=ceil(abs(imax-imin)+1);
% axis([imin imax 0-jmax 0-jmin]);
if nargin < 11 || nrows < 1
    nrows=ceil(abs(iimax-imin)+1)
end
if nargin < 12 || ncols < 1
    ncols=ceil(abs(jjmax-jjmin)+1)
end
%axis([jmin jmax 0-imax 0-imin]);
% 2014-JAN-07 SWAP ROWS AND COLUMNS
% iimin
% iimax 
% jjmin 
% jjmax
axis([iimin iimax 0-jjmax 0-jjmin]);

axis xy;
% make one patch as big as the image
% xtemp =     [jmin jmax+1 jmax+1 jmin  ];
% ytemp =   0-[imin imin   imax+1 imax+1];
% 2014-JAN-07 SWAP ROWS AND COLUMNS
xtemp =     [iimin iimax+1 iimax+1 iimin  ];
ytemp =   0-[jjmin jjmin   jjmax+1 jjmax+1];
ctemp = cmap(1,:);

if drawit == 1
    patch(xtemp,ytemp,ctemp,'LineStyle','none');
end
imout = ones(nrows,ncols);
% plot missing data in gray
%imout = floor(numel(cuts)/2)*ones(nrows,ncols);
%imout = nan(nrows,ncols);
i0 = floor(min(qi1));
j0 = floor(min(qj1));


nlevels = numel(cuts);

% slow, but works
for i=1:numel(qqq)
    % coordinates of patch corners
    %     xtemp =    [qj1(i) qj2(i)+1 qj2(i)+1 qj1(i)  ];
    %     ytemp = 0- [qi1(i) qi1(i)   qi2(i)+1 qi2(i)+1];
    % 2014-JAN-07 SWAP ROWS AND COLUMNS
    xtemp =    [qi1(i) qi2(i)+1 qi2(i)+1 qi1(i)  ];
    ytemp = 0- [qj1(i) qj1(i)   qj2(i)+1 qj2(i)+1];
    
    % indices
    ii=(qi1(i):qi2(i))-i0+1;
    jj=(qj1(i):qj2(i))-j0+1;
    
    % choose the color
    qtemp = qqq(i);
    ctemp = NaN;
    for j=1:nlevels-1
        if qtemp >= cuts(j) && qtemp < cuts(j+1)
            ctemp = cmap(j,:);
            % index
            %imout(ii,jj)= j;
            % 2014-JAN-07 SWAP ROWS AND COLUMNS
            imout(jj,ii)= j;
        end
    end
    
    % draw the patch
    if drawit == 1
        h=patch(xtemp,ytemp,ctemp,'LineStyle','none');
    end
end

% draw the patches
%h=patch(xtemp,ytemp,ctemp,'LineStyle','none');

if drawit == 1
    title(tstr);
    xlabel('column index');
    ylabel('row index');
    
    % label color bar with cutting values
    if drawcolorbar == 1
        
        for i=1:nlevels-1
            ii=find((isfinite(qqq)==1) & (qqq >= cuts(i)) & (qqq < cuts(i+1)));
            %fprintf(1,'%5d %+10.4e  %+10.4e %10d\n',i,cuts(i),cuts(i+1),numel(ii));
        end
        
        h=colorbar;
        set(h,'YTickMode'     ,'Manual');
        set(h,'YTickLabelMode','Manual');
        set(h,'Clipping'      ,'Off');
        
        YTick=get(h,'YTick');
        YTick = [0:8:nlevels]/nlevels;
        set(h,'YTick',YTick);
        
        for i=1:numel(YTick)-1
            YTickLabel{i} = sprintf('%10.4f',cuts(8*(i-1)+1));
            %fprintf(1,'%d %f %s\n',i,YTick(i),YTickLabel{i});
        end
        YTickLabel{numel(YTick)}= sprintf('%10.4f',nanmax(nanmax(qqq)));
        set(h,'YTickLabel',YTickLabel);
        ylabel(h,zstr);
    end
end

% replace nan values with a 1 to point to the lowest value in the color
% table
% inan = find(isnan(imout)==1);
% imout(inan)= 1;

return;
end
