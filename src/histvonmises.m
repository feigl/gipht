function histvonmises(data,imode,titlestring)
%function histvonmises(data,imode)
%
% Produces a histogram of the data using Von Mises Distribution 
% input angular data should be in cycles in the range [-1/2, +1/2]
%
% imode == 0 do NOT remove mean
% imode == 1 REMOVE mean [default]
% 
%
% 
% For example, try
% 
% x = 0.2*pi+0.5*pi*randn([10000 1]);
% th = angle(complex(cos(x),sin(x)));
% histvonmises(th/2/pi)
% 
%
% Kurt Feigl
% last modified 2008-FEB-11

if nargin < 2
   imode = 1;
end

% prune
% iok = find(isreal(data)==1);
% iok = intersect(iok,find(isfinite(data)==1));
% data = data(iok);
% iok = find(isreal(data)==1);
% iok = intersect(iok,find(isfinite(data)==1));
% data = data(iok);
% remove complex

% 2011-NOV-30 take real part
if isreal(data) ~= 1
    warning(sprintf('Found complex values in data. Max abs(imaginary) is %g\n',max(imag(data))));
    data = real(data);
end

% prune 
iok = find(isfinite(data)==1);
data = data(iok);

% number of data points
n=numel(data);

if n < 50
    warning(sprintf('n (%d) too small\n',n));
    return
end
% number of bins
%nbins = floor(n/100);
%nbins = nbins - mod(nbins,2); % want an even number of bins

%if nbins < 10
%   nbins = 10;
%end
%if nbins > 50
%   nbins = 50;
%end
nbins=32;

width = 1/nbins; % width of bins in cycles

% if data are on half interval [0, 1/2], then double
if min(data) >= 0
   dbl = 1;
else
   dbl = 0;
end
data = data * (dbl+1);

% bin
%[Ncounts, BinCenters] = hist(data, nbins)
for i=1:nbins
   BinEdges(i) = -0.5 + (i-1)*width;
end
BinEdges(nbins+1)=0.5;
for i=1:nbins
   BinCenters(i) = (BinEdges(i)+BinEdges(i+1))/2;
end

 
%test if data are distributed as von Mises
% [Svm, Scrit95] = vonmisesness(data*2*pi);
% if Svm < Scrit95
%    fprintf(1,'Data are distributed as Von Mises because S = %g is less than Scrit95 = %g\n',Svm,Scrit95);
% else
%    fprintf(1,'Data are NOT distributed as Von Mises because S = %g is greater than Scri05 = %g\n',Svm,Scrit95);
% end


s0 = circular_stddev(data*2*pi); % in radians

md = mean_direction(data*2*pi);      % in radians
%R = 1 - circular_variance(data*2*pi)   % mean resultant length
Rbar = mean_resultant_length(data*2*pi);   % mean resultant length

d0 = circular_mean_deviation(data*2*pi,imode);   

if imode == 1
   % % shift mean to origin
   data2 = angle(complex(cos(data*2*pi - md),sin(data*2*pi - md)))/2/pi;
   xlab = 'difference from mean direction (cycles)';
else
   % do not shift
   data2 = data;
   xlab  = 'direction (cycles)';
end

% % shift mean to origin and convert to radians
data3 = angle(complex(cos(data2*2*pi - md),sin(data2*2*pi - md)));
% % Test if data are distributed as Von Mises - Mardia and Jupp
% [Sm, teststring1] = vonmisesness(data3);
% % Test if data are distributed as Von Mises - Fisher
% [U2, teststring2] = vonmisesness2(data3);


Ncounts = histc(data2, BinEdges);
%figure; 
%hist(data, nbins);  hold on;
bar(BinEdges,Ncounts,'histc');  hold on;
colormap([0.5 0.5 0.5]);
%axis([-0.5 0.5 0 ceil(10*1.1*max(Ncounts))/10]);
axis([-0.5 0.5 0 Inf]);
%axis([-0.5 0.5 0 ceil(10*n/4)/10]);
%axis([-0.5 0.5 0 ceil(10*n/15)/10]);
%axis([-0.5 0.5 0 ceil(10*n/6)/10]);

% for large n, then sigmastat is normally distributed
Cbar = sum(cos(data*2*pi))/n;
sigmastat = Cbar*sqrt(2*n);

% concentration parameter in Von Mises Distribution
kappa = batschelet_inv(Rbar);

%dcdf(1) = von_mises_cdf((-0.5+width)*pi,  md, kappa)-von_mises_cdf(-0.5*pi,  md, kappa);
for i=1:nbins
   e1 = BinCenters(i) - width/2;
   e2 = BinCenters(i) + width/2;
   if imode == 0
      dcdf(i) = von_mises_cdf(e2*2*pi,  md, kappa) ...
              - von_mises_cdf(e1*2*pi,  md, kappa);  % argument in radians
   else
      dcdf(i) = von_mises_cdf(e2*2*pi,  0, kappa) ...
              - von_mises_cdf(e1*2*pi,  0, kappa);  % argument in radians
   end
   % %   dcdf(i) = von_mises_cdf(e2*2*pi+md,  0, kappa) ...
   %           - von_mises_cdf(e1*2*pi+md,  0, kappa);  % argument in radians
end
%dcdf(1) = dcdf(end);
%dcdf(1:4)

plot(BinCenters,n*dcdf,'r*-');

% 
% % sort into quantiles
% x=sort(data2);
% 
% % initialize
% q=zeros(n,1);
% xqm=zeros(n,1);
% 
% for i=1:length(x)	
%   q(i)=(i-0.5)/n;
%   %xqm(i)=phiinv(q(i));
%   %xqm(i) = icdf(dist,q(i),a,b,c);
%   %xqm(i) = von_mises_cdf_inv(q(i), md, kappa)/2/pi; % in cycles 
%   xqm(i) = von_mises_cdf_inv(q(i),  0, kappa)/2/pi; % in cycles 
% end
% %h=figure;
% %h=plot(xqm,x,'r+');
% plot(xqm,x,'r+');
% hold on
% axis square
% axis ([-0.5 0.5 -0.5 0.5]);
% set(gca,'XTick',[-0.5:0.1:0.5]);
% set(gca,'YTick',[-0.5:0.1:0.5]);
% 
% %
% % Plot the theoretical line for the given md and kappa
% 
% %   h=plot([xqm(1) xqm(n)],[m+phiinv(0.5/n)*s m+phiinv((n-0.5)/n)*s ]  ,'k--');
% 
% y1 = (von_mises_cdf_inv(0.5/n,     0, kappa))/2/pi;
% y2 = (von_mises_cdf_inv((n-0.5)/n, 0, kappa))/2/pi;
% 
% plot ([xqm(1) xqm(n)],[y1 y2],'k-');
% legend('observed','ideal Von Mises','Location','NorthWest')
% 
% plot([0 0],[-0.5 +0.5],'k--');
% plot([-0.5 +0.5],[0 0],'k--');


xlabel(xlab,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
%ylabel('Number of occurences','FontName','Helvetica','Fontsize',14,'FontWeight','bold');
ylabel('n','FontName','Helvetica','Fontsize',14,'FontWeight','bold');
set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');

% statstring=sprintf(' n=%d DBL=%1d mu=%.3fcy d0=%.3fcy k=%.3f R=%.4f Sm=%.2f U2=%.4f'...
%     ,n,dbl,md/pi/2,d0/pi/2,kappa,Rbar,Sm,U2);

% statstring = strcat(statstring,'\n',teststring1);
% statstring = strcat(statstring,'\n',teststring2);

% if exist('titlestring','var') == 1
%     text(0.0,1.09,titlestring,'units','normalized','Clipping','Off','FontName','Courier','Fontsize',7,'FontWeight','normal');
% end
% text(0.0,1.05,teststring1,'units','normalized','Clipping','Off','FontName','Courier','Fontsize',7,'FontWeight','normal');
% text(0.0,1.01,teststring2,'units','normalized','Clipping','Off','FontName','Courier','Fontsize',7,'FontWeight','normal');

% if nargin < 3
%     title(statstring,'FontName','Courier','Fontsize',10,'FontWeight','normal');
% else
%     title(strcat(titlestring, statstring),'FontName','Courier','Fontsize',10,'FontWeight','normal');
% end
% hold on; plot(0,0,'.w'); % draw white dot
% legend(titlestring,statstring,teststring,'Location','NorthEastOutside')


return


