function [Xlag,Ylag,Rlag,Acov,Corr] = analyze_covariance2(X,Y,V)
% analyze and plot covariance of a 2-D field
% inputs
%    X,Y coordinates in meters
%    V   field value in meters
% outputs
% 2021/09/18 Kurt Feigl
%
nf=0;

% find extrema
xmin=nanmin(colvec(X));
xmax=nanmax(colvec(X));
ymin=nanmin(colvec(Y));
ymax=nanmax(colvec(Y));

% prune
iok=find(isfinite(V));
iok=intersect(iok,find(isfinite(X)));
iok=intersect(iok,find(isfinite(Y)));
X=X(iok);
Y=Y(iok);
V=V(iok);


% generate extrapolating function
Fs=scatteredInterpolant(X,Y,V,'natural','nearest');

% make a rectangular array for autocovariance
% step size
% dx=nanmean(diff(unique(colvec(X))));
% dy=nanmean(diff(unique(colvec(Y))));
dx=(xmax-xmin)/500
dy=(ymax-ymin)/500
xvec=[xmin:dx:xmax];
yvec=[ymin:dy:ymax];
[Xgrd,Ygrd] = meshgrid(xvec,yvec);

% evaluate interpolating function at each node of grid
Darray=Fs(Xgrd,Ygrd);
whos Darray

% plot the interpolated field
figure;
imagesc(xvec,yvec,Darray);
hold on;
plot(X,Y,'ko');
axis xy;
axis equal;
axis tight;
colormap(jet);
colorbar;
xlabel('X');
ylabel('Y');
title('Interpolated field');

% autocovariance
Acov=xcorr2(Darray);
[nDcorr,mDcorr] = size(Acov);
whos Acov

% lag vectors
xlag=[-(mDcorr-1)/2:(mDcorr-1)/2]*dx;
ylag=[-(nDcorr-1)/2:(nDcorr-1)/2]*dy;


% Correlation coefficient
%Corr =  inv(sqrtm(diag(Acov))) * Acov * inv(sqrtm(diag(Acov)));
Corr = Acov ./ diag(diag(Acov)); 
% histogram of autocovariance
nf=nf+1;figure
histogram(colvec(Acov));
xlabel('autocovariance');
ylabel('count');

% map of autocovariance
figure
imagesc(xlag,ylag,Acov);
colormap(jet);
axis xy;
axis equal;
axis tight;
xlabel('lag in X');
ylabel('lag in Y');
colorbar;
title('Autocovariance');
%printpng(sprintf('%s_autocovariance.png', mfilename));

% Distance
[Xlag,Ylag] = meshgrid(xlag,ylag);
Rlag=hypot(Xlag,Ylag);
whos Xlab Ylag Rlag
nf=nf+1;figure;
imagesc(xlag,ylag,Rlag);
colormap(jet);
axis xy;
axis equal;
axis tight;
xlabel('lag in X');
ylabel('lag in Y');
colorbar;
title('Distance ');

% plot orrelation coefficient as a function of distance
nf=nf+1;figure;
%loglog(colvec(Rlag)/1.e3,1.e3*sqrt(colvec(Dcorr)),'.r');grid on;
%plot(colvec(Rlag)/1.e3,1.e3*sqrt(colvec(Dcorr)),'.r');
plot(colvec(Rlag),colvec(Corr),'.r');
xlabel('lag distance');
%ylabel('mm');
ylabel('Corr');
%title('sqrt(autocovariance)');
title('Correlation coefficient');
%printpng(sprintf('%s_AcovDist.png', mfilename));
return
end

