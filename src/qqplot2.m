function [h,phat,chi2gof_h,chi2gof_p,chi2gof_stats] = qqplot2(data,dist)
%function [h,phat,chi2gof_h,chi2gof_p,chi2gof_stats] = qqplot2(data,dist)
%
%
% Produces a Quantile Quantile plot of the data,
%     and, optionally, perform a chi-squared test of goodness of fit to the
%     distribution
%
% inputs
%      data
%      dist  == name of distribution
% outputs
%      h     == plot handle of QQ poot
%      phat  == Maximum Likelihood estimate of parameters in distribution
%      normr == a measure of goodness of fit
%      chi2gof_h, chi2gof_p, chi2gof_stats == optional output: see help chi2gof
%
%
% Roots: example B.8, pages 266-267
% Aster, Thurber  and Borcher
%
% modified by Kurt Feigl 2011-03-31
% revised 20140107 Kurt Feigl to calculate correlation coefficient properly
% 2020/05/28 changed name to qqplot2
%
% Calls stats tool box routines icdf and mle
% For a list of permissable values of 'dist', type:
%     help mle
%
% Examples:
% 
% aa=randn(1000,1); % generate a sample of 1000 data 
%                   from a normal distribution
%
% qqplot2(aa,'normal');
%
%    The Null Hypothesis cannot be rejected at the 5 percent significance level.
%    Interpretation: sample is compatible with a normal distribution.
%
% qqplot2(abs(aa),'Generalized Pareto');
%
% qqplot2(abs(aa),'exponential');
%


% initialize
h = NaN;
cor=NaN;
phat = NaN;
normr = NaN;
chi2gof_h = NaN;
chi2gof_p = NaN;
chi2gof_stats = NaN;

if license('checkout','Statistics_Toolbox') ~= 1
    warning(sprintf('Cannot find stats toolbox in %s\n',mfilename));
    return
end

if nargin < 2
    dist = 'norm';
end

n=numel(data);
data=reshape(data,n,1);
x=sort(data);
q=zeros(n,1);
xqm=zeros(n,1);


% get non-negative values
if   strcmpi(dist,'gamma')       == 1 ...
        || strcmpi(dist,'poisson')     == 1 ...
        || strcmpi(dist,'exponential')     == 1
    data = abs(data);
end
% rescale to positive values
if   strcmpi(dist,'lognormal')   == 1 ...
        || strcmpi(dist,'geometric')   == 1 ...
        || strcmpi(dist,'weibull')   == 1
    if min(data) < 0
        warning('Rescaling to positive values');
        data = 1.0 + (data - nanmin(data))/(nanmax(data)-nanmin(data));
    end
end

% estimated values of parameters
if   strcmpi(dist,'chisquare')   == 1
    %     MLE can also fit a custom distribution that you define using
    %     distribution functions, in one of three ways:
    %
    %     [...] = MLE(DATA,'pdf',PDF,'cdf',CDF,'start',START,...) returns MLEs
    %         for the parameters of the distribution defined by the probability
    %         density and cumulative distribution functions PDF and CDF.  PDF and CDF
    %         are function handles created using @.  They accept as inputs a vector
    %         of data and one or more individual distribution parameters, and return
    %         vectors of probability density values and cumulative probability
    %         values, respectively.  If the 'censoring' name/value pair is not
    %         present, you may omit the 'cdf' name/value pair.  MLE computes the
    %         estimates by numerically maximizing the distribution's log-likelihood,
    %         and START is a vector containing initial values for the parameters.
    %     phat = mle(data,'pdf',@(z)pdf(dist,z,phat(1)),);
    phat(1) = 1.0;
    m = 1;
else
    phat = mle(data,'distribution',dist);
    m = numel(phat);
end

% make bins for quantiles
for i=1:n
    % Centers of bins
    q(i)=(i-0.5)/n;
    % 2014-01-08 edges of bins, starting at zero
    %q(i)=(i-1)/n;
end
% if strcmpi(dist,'norm') == 1
%     xqm = phiinv(q);
% else
switch m
    case 1
        xqm = icdf(dist,q,phat(1));
        phat(2:3) = NaN;
    case 2
        xqm = icdf(dist,q,phat(1),phat(2));
        phat(3:3) = NaN;
    case 3
        xqm = icdf(dist,q,phat(1),phat(2),phat(3));
    otherwise
        m
        error('Number of parameters in distribution model too large');
end
% 2014-01-08 use left edge of bin
%xqm = xqm(1:n);
%end


% % draw best fitting line
[P,S,MU] = polyfit(xqm,x,1);
% normr = S.normr;
[Y,DELTA] = polyval(P,xqm,S,MU);


% perform goodness of fit test
fprintf(1,'Testing the null hypothesis that the data are random sample from a %s distribution\n',dist);
%     X is a random sample from a normal distribution)
%           H = CHI2GOF(X) performs a chi-square goodness-of-fit test that the data in
%     the vector X are a random sample from a normal distribution with mean and
%     variance estimated from X.  The result is H=0 if the null hypothesis (that
%     X is a random sample from a normal distribution) cannot be rejected at the
%     5% significance level, or H=1 if the null hypothesis can be rejected at
%     the 5% level.  CHI2GOF uses NBINS=10 bins, and compares the test statistic
%     to a chi-square distribution with NBINS-3 degrees of freedom, to take into
%     account that two parameters were estimated.

%         [H,P] = CHI2GOF(...) also returns the p-value P.  The P value is the
%     probability of observing the given result, or one more extreme, by
%     chance if the null hypothesis is true.  If there are not enough degrees
%     of freedom to carry out the test, P is NaN.

%        % Three equivalent ways to test against an unspecified normal
%        % distribution (i.e., with estimated parameters)
%        x = normrnd(50,5,100,1);
%        [h,p] = chi2gof(x)
%        [h,p] = chi2gof(x,'cdf',@(z)normcdf(z,mean(x),std(x)),'nparams',2)
%        [h,p] = chi2gof(x,'cdf',{@normcdf,mean(x),std(x)})

% significance level
%alpha = 1.0e-6;
alpha = 0.05; 
nbins = 20;
switch m
    case 1
        [chi2gof_h,chi2gof_p,chi2gof_stats] = chi2gof(data ...
            ,'cdf',@(z)cdf(dist,z,phat(1)),'nparams',m,'alpha',alpha,'nbins',nbins,'emin',5);
    case 2
        [chi2gof_h,chi2gof_p,chi2gof_stats] = chi2gof(data ...
            ,'cdf',@(z)cdf(dist,z,phat(1),phat(2)),'nparams',m,'alpha',alpha,'nbins',nbins,'emin',5);
    case 3
        [chi2gof_h,chi2gof_p,chi2gof_stats] = chi2gof(data ...
            ,'cdf',@(z)cdf(dist,z,phat(1),phat(2),phat(3)),'nparams',m,'alpha',alpha,'nbins',nbins,'emin',5);
    otherwise
        error(sprintf('Wrong number of parameters %d in distribution %s',m,dist));
end

if chi2gof_h == 0
    fprintf(1,'\nThe Null Hypothesis cannot be rejected at the %.1f percent significance level.\n',alpha*100);
    fprintf(1,'Interpretation: sample is compatible with a %s distribution.\n',dist);
else
    fprintf(1,'\nThe Null Hypothesis can be rejected at the %.1f percent significance level.\n',alpha*100);
end


%% now make histogram
% st = struct with fields:
%     chi2stat: 2.5550
%           df: 3
%        edges: [-0.5000 0.5000 1.5000 2.5000 3.5000 5.5000]
%            O: [6 16 10 12 6]
%            E: [7.0429 13.8041 13.5280 8.8383 6.0284]

% The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level. The vector E
% contains the expected counts for each bin under the null hypothesis, and O contains the observed counts for each bin.
binCenters = chi2gof_stats.edges(1:end-1) + diff(chi2gof_stats.edges/2.); % centers of bins
binEdges = chi2gof_stats.edges; % edges of bins, but last one is too wide
nbins = numel(binCenters);
valsObs = chi2gof_stats.O;
valsExp = chi2gof_stats.E;
ncounts  = numel(valsObs);

% calculate theoretical values from distribution
switch m
    case 1
        evals2 = n*cdf(dist,binEdges,phat(1));        
    case 2
        evals2 = n*cdf(dist,binEdges,phat(1),phat(2));       
    case 3
        evals2 = n*cdf(dist,binEdges,phat(1),phat(2),phat(3));
    otherwise
        error('Unknown case','m = %d',m);
end
% difference cdf
evals2 = diff(evals2);

% 
% fprintf(1,'  i %12s %12s %12s %12s %12s %12s\n','zvalsL','zvalsR','cvals','ovals','evals1','evals2');
% for i=1:nbins
%     fprintf(1,'%3d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',i,zvals(i),zvals(i+1),cvals(i),ovals(i),evals1(i),evals2(i));
% end

titlestr = sprintf('%s distribution H = %d (n = %d) P = %g phat %g %g %g '...
    ,dist,chi2gof_h,numel(data),chi2gof_p,phat(1:m));

%% draw QQ-plot
figure; hold on;
set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
plot(xqm,x,'ro');

%axis([min([x;xqm]) max([x;xqm]) min([x;xqm]) max([x;xqm])]);
%axis([0 1 0 1]);
axis square

plot(xqm,Y,'k-');
plot(xqm,Y-DELTA,'k--');
plot(xqm,Y+DELTA,'k--');
% plot([0 1],[0 1],'b:');

xlabel(sprintf('expected values'));
ylabel(sprintf('observed values'));

title(strcat('QQ plot for ',titlestr));
h(1) = gcf;

%% draw histogram
figure;hold on;
% plot using edges, but the right-most bin is too wide
for i=1:numel(binEdges)-1
    fill([binEdges(i) binEdges(i+1) binEdges(i+1) binEdges(i)],[0 0 valsObs(i) valsObs(i)],'y');
end

% plot using the built-in routine, with an odd number of bins
if mod(nbins,2) == 0
    nbins2 = nbins+1;
else
    nbins2 = nbins;
end
histogram(data,nbins2);

% plot expected values at center 
plot(binCenters,valsExp,'ro-','MarkerSize',6,'LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k');
set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
title(sprintf('%s\nBlue:Observed, Yellow:Chi2gof, Red:expected',titlestr));
xlabel('Value');
ylabel('Number of occurrences');
h(2) = gcf;


%% draw another QQ plot using bins from chi2gof 
figure; hold on;
[dummy,isort] = sort(valsExp);
plot(valsExp(isort)/n,valsObs(isort)/n,'ro-');
axis([0 1 0 1]);
axis square;
hold on;
plot([0 1],[0 1],'b:');
set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
title(titlestr);
xlabel('Expected frequency');
ylabel('Observed frequency');
h(3) = gcf;

%% draw residual 
figure; hold on;
plot(valsExp(isort)/n,(valsObs(isort)-valsExp(isort))/n,'ro-');
%axis([0 1 -Inf +Inf]);
axis tight
hold on;
plot([0 nanmax(valsExp(isort))/n],[0 0],'b:');
title(titlestr);
ylabel('Observed - expected frequency');
xlabel('Expected frequency');

h(4) = gcf;

return

% function x=phiinv(z)
% %function x=phiinv(z)
% %
% %  Calculates the inverse normal distribution.
% %
% %   x=phiinv(z)
% %
% %  gives z=int((1/sqrt(2*pi))*exp(-t^2/2),t=-infinity..x)
% %
% if (z >= 0.5),
%     x=sqrt(2)*erfinv((z-0.5)/.5);
% else
%     x=-phiinv(1-z);
% end
% return
% 



