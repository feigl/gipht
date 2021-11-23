function [Xlag,Ylag,Rlag,Acov,Corr1] = analyze_covariance(X,Y,V)
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

% calculate distance from center
%R=hypot(X-nanmean(X),Y-nanmean(Y));
%calculate distance 
R=hypot(diff(X),diff(Y));
% sort by distance
[R,isort]=sort(R);
whos R
min(R)
V=V(isort);
ndata=numel(V)

figure
plot(R,'g-');
xlabel('index');
ylabel('distance R');


% autocovariance
%[Acov,lags1]=xcov(V);
%[Acov,lags1]=xcov(V,'unbiased'); % 
% [nDcorr,mDcorr] = size(Acov);
% whos Acov
[Corr1,lags1]=xcov(V,'coeff'); % 


% Correlation coefficient
% normalize by autocovariance at zero lag
%Corr1 = Acov / Acov(ndata);
whos Corr1

% Use Matlab function
% no scaling
%[Corr2,lags2]  = xcorr(V,'none'); 
% normalize the sequence so that the auto-correlations at zero lag are identically 1.0.
[Corr2,lags2]  = xcorr(V,'coeff'); 
whos Corr2

% Calculate 
%Lohman, R. B., and M. Simons (2005), Some thoughts on the use of InSAR data to constrain models of surface deformation: Noise structure and data downsampling, Geochemistry, Geophysics, Geosystems, 6, n/a-n/a. 


% % histogram of autocovariance
% nf=nf+1;figure
% histogram(colvec(Acov));
% xlabel('autocovariance');
% ylabel('count');

% histogram of correlation
nf=nf+1;figure
histogram(colvec(Corr1));
xlabel('Correlation coefficient');
ylabel('count');

% Correlation as a function of distance
nf=nf+1;figure; hold on;
% plot(lags1(ndata:end),Corr1(ndata:end),'r.');
% plot(abs(lags1(1:ndata)),Corr1(1:ndata),'ro');
plot(lags2(ndata:end),Corr2(ndata:end),'b.');
plot(abs(lags2(1:ndata)),Corr2(1:ndata),'bo');
%legend('Corr1 right','Corr1 left','Corr2 right', 'Corr2 left');
grid on;
xlabel('lag');
ylabel('correlation');
title('Correlation coefficient');

nf=nf+1;figure; hold on;
% plot(R(2:end),Corr1(ndata+1:end),'r.');
% plot(R(2:end),Corr1(ndata-1:-1:1),'ro');
plot(R(2:end),Corr2(ndata+1:end),'b.');
plot(R(2:end),Corr2(ndata-1:-1:1),'bo');
%legend('Corr1 right','Corr1 left','Corr2 right', 'Corr2 left');
grid on;
xlabel('distance [m]');
ylabel('correlation');
title('Correlation coefficient');

return
end

