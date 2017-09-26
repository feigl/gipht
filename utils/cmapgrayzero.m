function cmap = cmapgrayzero
% make a color map with zero set to gray
cmap=colormap('jet');
[nlevels,ndum] = size(cmap)
%cmap=zeros(nlevels+1,3);
% cmap=zeros(nlevels,3);
% for i=1:nlevels-1
%     cmap(i+1,:)=tmap(i,:);
% end
i1 = floor(nlevels/2)-1;
i2 = ceil(nlevels/2)+1;
for i=i1:i2
    cmap(i,1:3) = [0.5 0.5 0.5]; % gray
end
% for i=1:nlevels
%     fprintf(1,'%3d %10.4f %10.4f %10.4f\n',i,cmap(i,1),cmap(i,2),cmap(i,3));
% end
%colormap(cmap)
return;
end


