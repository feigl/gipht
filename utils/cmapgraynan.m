function cmap=cmapgraynan
% make a color map with NaN set to gray
tmap= colormap('jet');
[nlevels,ndum] = size(tmap);
cmap=zeros(nlevels+1,3);
for i=1:nlevels
    cmap(i+1,:)=tmap(i,:);
end
cmap(1,:) = [0.5 0.5 0.5]; % gray in R,G,B from 0 to 1
colormap(gca,cmap);
return;
