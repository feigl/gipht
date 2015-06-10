function crit69 = test_generalized_paretos(sample1,mean1)
%function P = testGPmeans(gpmean1, n1, gpmean2, n2)
% Return critical value of mean for sample2 drawn from
% a Generalized Pareto (GP) distribution
% such that sample2 is significantly different from sample1
%
% Kurt Feigl 2011-OCT-11
%
% INPUTS:
%    gpmean1, gpmean2 GP of each sample, as estimated by gpfit
%    n1, n2 number of data in each sample
% OUTPUT:
%    eta == test statistic given by 7.3.23
%           it is normally distributed as N(0,1)
%           zero mean and unit variance
%           The critical region is both tails.
%           We accept Null Hypothesis H0 with 95 % confidence if
%               abs(eta) > 1.96
%
%
%
% Test Case for n1 = n2 = 10
%    P=testkappas3(0.2,10,0.1,10)
%    eta = icdf('normal',P,0,1)
% so eta = eta1 QED

% Test assumptions
if nargin ~= 2
    error('Need 2 arguments');
end

n = numel(sample1);
if n < 10
    warning(sprintf('Sample size %d is too small',n));
end


%     help gpfit
%     GPFIT Parameter estimates and confidence intervals for generalized Pareto data.
%     PARMHAT = GPFIT(X) returns maximum likelihood estimates of the parameters
%     of the two-parameter generalized Pareto (GP) distribution given the data
%     in X.  PARMHAT(1) is the tail index (shape) parameter, K and PARMHAT(2) is
%     the scale parameter, SIGMA.  GPFIT does not fit a threshold (location)
%     parameter.
%

% estimate parameters for Generalized Pareto Distribution
[phat psig]  = gpfit(sample1,0.05);
if numel(isfinite(phat)) ~= 2 || numel(isfinite(psig)) ~= 4
    warning('Assumptions not valid');
    phat
    psig
    crit69 = NaN;
    return;
end

% Uncorrelated noise in radians by drawing from GP distribution
% GPRND Random arrays from the generalized Pareto distribution.
%     R = GPRND(K,SIGMA,THETA) returns an array of random numbers chosen from the
%     generalized Pareto (GP) distribution with tail index (shape) parameter K,
%     scale parameter SIGMA, and threshold (location) parameter THETA.  The size
%     of R is the common size of K, SIGMA, and THETA if all are arrays.  If any
%     parameter is a scalar, the size of R is the size of the other parameters.
%
%     R = GPRND(K,SIGMA,THETA,M,N,...) or R = GPRND(K,SIGMA,[M,N,...]) returns
%     an M-by-N-by-... array.
theta = 0.0;
GP = gpmv(phat);

kount = 0;
h = 0;
mean2 = GP.m;

% brute force grid search until significantly different
while h == 0 && kount < 100
    kount = kount + 1;
    mean2 = 1.01 * mean2;
    phat(1) = (GP.a/mean2) - 1.0;
    sample2 = gprnd(phat(1),phat(2),theta,size(sample1)); % column vector
    
    %  KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
    %     H = KSTEST2(X1,X2,ALPHA,TYPE) performs a Kolmogorov-Smirnov (K-S) test
    %     to determine if independent random samples, X1 and X2, are drawn from
    %     the same underlying continuous population. ALPHA and TYPE are optional
    %     scalar inputs: ALPHA is the desired significance level (default = 0.05);
    %     TYPE indicates the type of test (default = 'unequal'). H indicates the
    %     result of the hypothesis test:
    %        H = 0 => Do not reject the null hypothesis at significance level ALPHA.
    %        H = 1 => Reject the null hypothesis at significance level ALPHA.
    %    Let S1(x) and S2(x) be the empirical distribution functions from the
    %     sample vectors X1 and X2, respectively, and F1(x) and F2(x) be the
    %     corresponding true (but unknown) population CDFs. The two-sample K-S
    %     test tests the null hypothesis that F1(x) = F2(x) for all x, against the
    %     alternative specified by TYPE:
    %         'unequal' -- "F1(x) not equal to F2(x)" (two-sided test)
    %         'larger'  -- "F1(x) > F2(x)" (one-sided test)
    %         'smaller' -- "F1(x) < F2(x)" (one-sided test)
    
    alpha = 0.31;
    type = 'larger';
    h = kstest2(sample1,sample2,alpha,type);
    fprintf(1,'kount = %3d mean2 = %20.8f h = %d\n',kount,mean2,h)
end
if h == 1
    crit69 = mean2;
else
    crit69 = NaN;
    return
end
end

function GP = gpmv(phat)
% given esimated parameters for generalized Pareto Distribution,
% return mean and variance
% Reference:
% Choulakian, V. and M.A. Stephens (2001),
%     Goodness-of-Fit Tests for the Generalized Pareto Distribution
%     Technometrics, 43, pp 478-485
%
GP.k = phat(1); % PARMHAT(1) is the tail index (shape) parameter, K
GP.a = phat(2); % PARMHAT(2) is the scale parameter, SIGMA.

if abs(1.0 + GP.k) > 0.0 && abs(1 + 2.0 * GP.k) > 0.0
    GP.m = GP.a /(1.0 + GP.k);  % mean mu
    GP.v = GP.a^2/(((1.0 + GP.k)^2) * (1 + 2.0 * GP.k)); % variance
else
    GP.m = NaN;
    GP.v = NaN;
end
return
end



