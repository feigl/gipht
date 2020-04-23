%https://www.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc
close all
clear all
x=[-5:+5];
y=[-10:+10];
[X,Y] = meshgrid(x,y);
data = X .* Y;
%data = rand(10);
data(X > 4) = NaN;
[nr,nc] = size(data);
subplot(3,1,1);
imagesc(data);
colorbar;
subplot(3,1,2);
alpha=ones(size(data));
alpha(isnan(data))=0;
imagesc(x,y,data,'AlphaData',alpha);
%alpha=true(size(data));
% alpha=isnan(data);
% imagesc(x,y,data,'AlphaData',alpha,'AlphaDataMapping','direct');
colorbar;
subplot(3,1,3);
pcolor(X,Y,data);
colorbar;
shading flat;
colormap(jet);
%set(gca, 'ydir', 'reverse');

