function Rbar1 = testkappas_inv3(P, n1, Rbar2, n2)
%
% given P, n1, Rbar2, and n2,
%
% find corresponding Rbar1
%
%Test equality of concentration parameters
% Mardia and Jupp, Directional Statistics, Wiley, 1999
% chap. 7.3.2 Tests of Equality of Concentration Parameters, pages 133-134
% Null hypothesis is that concentration parameters are equal kappa1 = kappa2
% INPUTS:
%    Rbar1, Rbar2 mean resultant length of each sample
%    n1, n2 number of data in each sample
% OUTPUT:
%    eta == test statistic given by 7.3.23
%           it is normally distributed as N(0,1) 
%           zero mean and unit variance
%           The critical region is both tails.
%           We accept Null Hypothesis H0 with 95 % confidence if 
%               abs(eta) = 1.96
% 
%
% Kurt Feigl 2008-FEB-10
%

% Test assumptions
if nargin ~= 4
    error('Need 4 arguments')
end
if n1 < 5 || n2 < 5
    warning(sprintf('Sample sizes %d %d too small.\n',n1,n2));
end

if isfinite(Rbar2) ~= 1 || isreal(Rbar2) ~= 1
    Rbar1 = NaN;
    warning(sprintf('Mean resultant length Rbar2 %10.5f must be finite and real.\n',Rbar2));
    return;
end

% if Rbar2 > 0.70
%     Rbar1 = NaN;
%     warning(sprintf('Mean resultant length Rbar2 %10.5f too large for Case I or Case II\n',Rbar2));
%     return;
% end

% eta is supposed to be normally distributed,
%  so let's restrict to +/- 10 sigma
% if eta < -10 | eta > 10
%     Rbar1 = NaN;
%     warning(sprintf('Normally distributed eta %e is extreme!\n',eta));
%     return;
% end

% if Rbar2 > 0.99
%     guess = 2.e-6;
% elseif Rbar2 > 0.94
%     guess = 0.00002;
% elseif Rbar2 > 0.9
% %if Rbar2 > 0.9
%     guess = 0.0002;
% elseif Rbar2 > 0.7
%     guess = 0.02;
% elseif Rbar2 > 0.45
%     guess = 0.15;
% else
%     guess = 0.2;
% end
%guess = 0.42;
%Rbar1 = fzero(@(x)testkappas1(x, n1, Rbar2, n2)-eta,guess);

% 2012-JAN-03
if Rbar2 > 0 && Rbar2 < 1.00
    guess = Rbar2;
else
    guess = 0.5;
end
Rbar1 = fzero(@(x)testkappas3(x, n1, Rbar2, n2)-P,guess);
return

