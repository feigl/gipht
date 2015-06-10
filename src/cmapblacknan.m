function cmap=cmapblacknan
% make a color map with NaN set to black
tmap= colormap('jet');
[nlevels,ndum] = size(tmap);
cmap=zeros(nlevels+1,3);
for i=1:nlevels
    cmap(i+1,:)=tmap(i,:);
end
cmap(1,:) = [0 0 0]; % black
colormap(gca,cmap);
return;
