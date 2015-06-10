function cmap=cmapblackzero(ival)
% make a color map with zero set to black
if nargin > 0
    tmap= colormap('jet');
    [nlevels,ndum] = size(tmap);
    %cmap=zeros(nlevels+1,3);
    cmap=zeros(nlevels,3);
    for i=1:nlevels-1
        cmap(i+1,:)=tmap(i,:);
    end
    cmap(1,:)       = [0 0 0]; % black
    %cmap(nlevels,:) = [1 1 1]; % white
    %cmap(floor(nlevels/2),:) = [0.5 0.5 0.5]; % gray
else
    cmap=colormap('jet');
    cmap(33,:)=[0 0 0];
end
%colormap(gca,cmap);
return;
