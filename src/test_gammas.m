function crit69 = test_gammas(sample1)
%function crit69 = test_gammas(sample1)
% Return critical value of mean for sample2 drawn from
% a Gamma distribution
% such that sample2 is significantly different from sample1
%
% When a is large, the gamma distribution closely approximates a normal
% distribution with the advantage that the gamma distribution has density
% only for positive real numbers.

% Kurt Feigl 2014-JAN-08
%
% INPUTS:

% Test assumptions
if nargin ~= 1
    error('Need 1 argument');
end

n = numel(sample1);
if n < 10
    warning(sprintf('Sample size %d is too small',n));
end



% estimate parameters for gammma Distribution
% %     help gamfit
% gamfit Parameter estimates for gamma distributed data.
%     PARMHAT = gamfit(X) returns maximum likelihood estimates of the
%     parameters of a gamma distribution fit to the data in X.  PARMHAT(1)
%     and PARMHAT(2) are estimates of the shape and scale parameters A and B,
%     respectively.
% The mean of the gamma distribution with parameters a and b is ab. The variance is ab^2.

[phat psig]  = gamfit(sample1,0.05);
if numel(isfinite(phat)) ~= 2 || numel(isfinite(psig)) ~= 4
    warning('Assumptions not valid');
    phat
    psig
    crit69 = NaN;
    return;
else
    gamma1A = phat(1);
    gamma1B = phat(2);
end

% Uncorrelated noise in radians by drawing from GP distribution
%  gamrnd Random arrays from gamma distribution.
%     R = gamrnd(A,B) returns an array of random numbers chosen from the
%     gamma distribution with shape parameter A and scale parameter B.  The
%     size of R is the common size of A and B if both are arrays.  If
%     either parameter is a scalar, the size of R is the size of the other
%     parameter.
%  
%     R = gamrnd(A,B,M,N,...) or R = gamrnd(A,B,[M,N,...]) returns an
%     M-by-N-by-... array.
%  
%     Some references refer to the gamma distribution with a single
%     parameter.  This corresponds to gamrnd with B = 1.
%  

% brute force grid search until significantly different
%var2 = gamma1A * gamma1B^2;

kount = 0;
h = 0;
mean2 = nanmean(sample1);
mean2 = gamma1A * gamma1B;

while h == 0 && kount < 100
    kount = kount + 1;
%     var2 = 1.005 * var2;
%         gamma2A = var2/gamma1B^2;
%         gamma2B = sqrt(var2/gamma1A);
    mean2 = (1.0 + 5.0e-3) * mean2;
    gamma2A = mean2/gamma1B;
    gamma2B = mean2/gamma1A;
       
    sample2 = gamrnd(gamma2A,gamma2B,size(sample1)); % column vector
    
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
    %fprintf(1,'kount = %3d var2 = %20.8f h = %d\n',kount,var2,h)
    fprintf(1,'kount = %3d mean2 = %20.8f h = %d\n',kount,mean2,h)
end
if h == 1
    % return L1 norm
    %crit69 = nanmean(sample2);
    crit69 = mean2;
else
    crit69 = NaN;
    return
end
end


