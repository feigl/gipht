%function qqplotvonmises(data,titlestring)
%
% Produces a QQ plot of the data using Von Mises Distribution 
% input angular data should be in cycles in the range [-1/2, +1/2]
%
% Kurt Feigl
% last modified 2008-MAR-30

function qqplotvonmises(data,titlestring)

% number of data points
n=length(data);

if max(max(data)) - min(min(data)) > 1.0
    warning('Difference between extreme values exceeds 1.0 cycles');
end

% if data are on half interval [0, 1/2], then double
if min(data) >= 0
   dbl = 1;
else
   dbl = 0;
end
data = data * (dbl+1);

s0 = circular_stddev(data*2*pi); % in radians 

md = mean_direction(data*2*pi);      % in radians
%R = 1 - circular_variance(data*2*pi)   % mean resultant length
Rbar = mean_resultant_length(data*2*pi);   % mean resultant length

% for large n, then sigmastat is normally distributed
Cbar = sum(cos(data*2*pi))/n;
sigmastat = Cbar*sqrt(2*n);

% concentration parameter in Von Mises Distribution
kappa = batschelet_inv(Rbar);

d0 = circular_mean_deviation(data*2*pi,1);   

% shift mean to origin, keep in cycles
data2 = angle(complex(cos(data*2*pi - md),sin(data*2*pi - md)))/2/pi;

% % shift mean to origin and convert to radians
data3 = angle(complex(cos(data2*2*pi - md),sin(data2*2*pi - md)));
% Test if data are distributed as Von Mises - Mardia and Jupp
[Sm, teststring1] = vonmisesness(data3);
% Test if data are distributed as Von Mises - Fisher
[U2, teststring2] = vonmisesness2(data3);



% sort into quantiles
x=sort(data2);

% initialize
q=zeros(n,1);
xqm=zeros(n,1);


for i=1:length(x)	
  q(i)=(i-0.5)/n;
  %xqm(i)=phiinv(q(i));
  %xqm(i) = icdf(dist,q(i),a,b,c);
  %xqm(i) = von_mises_cdf_inv(q(i), md, kappa)/2/pi; % in cycles 
  xqm(i) = von_mises_cdf_inv(q(i),  0, kappa)/2/pi; % in cycles 
end



plot(xqm,x,'r+');
hold on
axis square
axis ([-0.5 0.5 -0.5 0.5]);
set(gca,'XTick',[-0.5:0.1:0.5]);
set(gca,'YTick',[-0.5:0.1:0.5]);
set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');

 
% Plot the theoretical line for the given md and kappa

%   h=plot([xqm(1) xqm(n)],[m+phiinv(0.5/n)*s m+phiinv((n-0.5)/n)*s ]  ,'k--');

y1 = (von_mises_cdf_inv(0.5/n,     0, kappa))/2/pi;
y2 = (von_mises_cdf_inv((n-0.5)/n, 0, kappa))/2/pi;

plot ([xqm(1) xqm(n)],[y1 y2],'k-','LineWidth',2);
h=legend('observed','ideal','Location','NorthWest');
set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');

plot([0 0],[-0.5 +0.5],'k--');
plot([-0.5 +0.5],[0 0],'k--');

xlabel('von Mises quantiles','FontName','Helvetica','Fontsize',14,'FontWeight','bold');
ylabel('sample quantiles'   ,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
% title(sprintf('DBL=%1d Mean.Dir = %.3f cy Cir.Std.Dev = %.3f cy Kappa = %.3f Rbar = %8.2f sigmastat = %.2f',dbl,md/pi/2,s0/pi/2,kappa,Rbar,sigmastat));
statstring=sprintf(' DBL=%1d mu=%.3fcy d0=%.3fcy k=%.3f R=%.4f Sm=%.2f U2=%.4f'...
   ,dbl,md/pi/2,d0/pi/2,kappa,Rbar,Sm,U2);
if nargin < 3
   title(statstring);
 else
    title(strcat(titlestring, statstring));
%    text(-0.48,1.15*max(Ncounts),statstr);
end

return


