% Case 1
n = 50;
x = 1:n;
y = 1:n;
data = peaks(n).^4;
contour_levels = [0.5, 10, 100, 500, 1000 ];
contour_levels2 = [0,contour_levels,5000];

% % Case 2
% x=1:0.5:5;
% y=1:0.5:10;
% [x,y] = meshgrid(x,y);
% data = x.^2+exp(y);
% data(end-5:end,end-3:end) = nan;
% data(end-2:end,end-2:end) = inf;
% contour_levels = [5,20,100,1000,10000];
% contour_levels2 = [0,contour_levels,50000];


figure('position',[1,1,1150,800])
subplot(331)
contourf(x,y,data)
colorbar,colormap(jet)
title('original linear cmap')
subplot(332)
contourfnu(x,y,data,contour_levels)
title('imagesc default')
subplot(333)
contourfnu(x,y,data,contour_levels,[],[],[],'contourf')
title('contourf default')
subplot(334)
contourfnu(x,y,data,contour_levels,[],[],false)
title('imagesc overticklabel=false')
subplot(335)
cmap = jet;
cmap(1,:) = [1,1,1]; %white
contourfnu(x,y,data,contour_levels,cmap,[],false)
title('imagesc FirstColorWhite')
subplot(336)
contourfnu(x,y,data,contour_levels2)
title('imagesc AnotherContourLevels')
subplot(337)
contourfnu(x,y,data,contour_levels,[],[],[],[],2)
title('imagesc interp')
subplot(338)
contourfnu(x,y,data,[-inf,contour_levels2,inf],[],[],false)
title('imagesc InfInLevels')
subplot(339)
contourfnu(x,y,data,contour_levels2,[],'southoutside')
title('imagesc PosCbar=''southoutside''')
