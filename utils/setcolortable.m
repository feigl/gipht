function setcolortable(ival)
% make a color map 
gray = [0.5 0.5 0.5]; % gray
black = [0. 0. 0.]; % black
white = [1. 1. 1.]; % white

if nargin < 1
    ival = 0;
end

switch ival
    case 0
        rgb = black;
    case 1
        rgb = white;
    otherwise
        rgb = gray;
end
cmap=colormap('jet');
[nlevels,ndum] = size(cmap);
%cmap=zeros(nlevels+1,3);
% cmap=zeros(nlevels,3);
% for i=1:nlevels-1
%     cmap(i+1,:)=tmap(i,:);
% end
i1 = floor(nlevels/2)-1
i2 = ceil(nlevels/2)+1
for i=i1:i2
    cmap(i,1:3) = rgb;
end
% cmap(1,1:3) = white;
% cmap(nlevels,1:3) = black;
% for i=1:nlevels
%     fprintf(1,'%3d %10.4f %10.4f %10.4f\n',i,cmap(i,1),cmap(i,2),cmap(i,3));
% end
colormap(cmap)
return;
