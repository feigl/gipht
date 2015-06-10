function [S, teststring] = vonmisesness(theta1)
%function [S, teststring] = vonmisesness(theta)
% Test if directional data theta in radians follow Von Mises distribution
%
% Mardia and Jupp, Directional Statistics, Wiley, 1999
% chap. 7.5 Testing Von Misesness, pages 142-143
% (As printed in this textbook, the expressions contain two typographic
% errors, as noted below.)
%
% The null hypothesis is that theta is distributed as Von Mises
%
% This test was proposed by:
%
% Cox, D. R., Contribution to discussion of Mardia (1975a)
%  Journal of the Royal Statistical Society, Series B, 1975
% 
% Further details may be found in
% Barndorff-Nielsen,O.E.; Cox,D.R. 1989
% Asymptotic Techniques for Use in Statistics
% Chapman and Hall, London 
%
% This code is based on the expressions in 
% Edgeworth and Saddle-Point Approximations with Statistical Applications 
% O. Barndorff-Nielsen; D. R. Cox 
% Journal of the Royal Statistical Society. Series B (Methodological), Vol. 41, No. 3. (1979), pp. 
% 279-312. 
%
% 
% Example 1
%
% x = 0.2*pi*randn([1000 1]);
% th = angle(complex(cos(x),sin(x)));
% String  = vonmisesness(th)
% 
% Example 2
%
% th = randraw('vonmises', [pi/2, 3], [1 1e5]);
% String  = vonmisesness(th)
% 
% 
% Kurt Feigl 2008-MAR-17
% 2009-JAN-12 prune NaN
%
%format long
S = NaN;
teststring = 'undefined';

% prune NaN
iok = find(isfinite(theta1)==1);
theta = theta1(iok);
% prune
if isreal(theta1) == 0
    warning('data must be real');
    return
end



n = numel(theta);
mu = mean_direction(theta);     % in radians
Rbar = mean_resultant_length(theta);   % mean resultant length

% concentration parameter in Von Mises Distribution
k = batschelet_inv(Rbar);

%sc = sum(cos(2*(theta-mu)))/n -  I1(k)/I0(k) % B-N and Cox 1979
%sc = sum(cos(2*(theta-mu)))/n -  I2(k)/I0(k); % Peter Jupp, personal communiction 2008-JAN-22
 sc = sum(cos(2*(theta-mu))) -  n*I2(k)/I0(k); % Peter Jupp, personal communiction 2008-FEB-21

%ss = sum(sin(2*(theta-mu)))/n; 
 ss = sum(sin(2*(theta-mu))); % Peter Jupp, personal communiction 2008-FEB-21

vc1 = ((I0(k))^2 + I0(k)*I4(k) - 2*((I2(k))^2))/(2*(I0(k))^2);

vc2 = ( I0(k)*I3(k) + I0(k)*I1(k) - 2*I1(k)*I2(k) )^2 ...
        / ...
       ( 2*(I0(k))^2 * ( (I0(k))^2 + I0(k)*I2(k) - 2 * (I1(k))^2 ));  % B-N and Cox 1979
%      ( 2*(I0(k))^2 * ( (I0(k))^2 + I0(k)*I1(k) - 2 * (I1(k))^2 ))  % Mardia and Jupp
  
vc = vc1 - vc2;

if vc < 0
    fprintf('WARNING: vc is negative %g\n',vc);
end

vs1 = I0(k) - I4(k);
vs2 = I0(k) - I2(k);
vs3 = (I1(k) - I3(k))^2;   % B-N and Cox 1979
%vs3 = (I0(k) - I3(k))^2;   % Mardia and Jupp
vs4 = 2*I0(k) * (I0(k) - I2(k)); % Peter Jupp, personal communiction 2008-JAN-22
vs = (vs1*vs2 - vs3) / vs4;
if vs < 0
    fprintf('WARNING: vs is negative %g\n',vs);
end

% Test statistic: Reject Null Hypothesis for large values of S
%S = sc^2/vc + ss^2/vs
%S = sc^2/vc/n/n + ss^2/vs/n/n % Chi-SQUARED
S1 = sc^2/vc/n + ss^2/vs/n; % Peter Jupp, personal communiction 2008-FEB-21                           
% S2 = (sc/vc/sqrt(n))^2 + (ss/vs/sqrt(n))^2 % Cox 1975 says "sum of squares, the latter distributed as Chi-2"
% S3 = (sc/vc/n)^2 + (ss/vs/n)^2 % Cox 1975 says "sum of squares, the latter distributed as Chi-2"
S = S1; % Finally works! 2008-MAR-19


% Under the null hypothesis, the large-sample asymptotic distribution of S
% is Chi-2 with 2 degrees of freedom
% The ICDF function is included in the Stats Toolbox
if length(which('icdf')) > 0
    Scrit95 = icdf('chi2',0.95,2);
    %Scrit05 = icdf('chi2',0.05,2);
else
    %Scrit95 = NaN;
    %warning('Cannot find STATISTICS TOOLBOX. Setting Scrit95 to NaN');
    Scrit95 = icdf0('chi2',0.95,2);
end

% The CDF function is included in the Stats Toolbox
if length(which('cdf')) > 0
    P = 1-cdf('chi2',S,2);
else
    P = 1-cdf0('chi2',S,2);
    %     P = NaN;
    %     warning('Cannot find STATISTICS TOOLBOX. Setting P to NaN');
end


if S < Scrit95
   teststring = sprintf('Data ARE compatible with Von Mises distribution because S = %g is less than critical Chi^2(2,0.05) = %g . (P = %8.2e)\n',S,Scrit95,P);
else
   teststring = sprintf('Data are NOT compatible with Von Mises distribution because S = %g is greater than critical Chi^2(2,0.05)  = %g . (P = %8.2e)\n',S,Scrit95,P);
end

if nargout < 2
   fprintf(1,'%s\n',teststring);
end

return


modified bessel functions of first kind order n
function i0 = I0(x)
i0=besseli(0,x);
return
function i1 = I1(x)
i1=besseli(1,x);
return
function i2 = I2(x)
i2=besseli(2,x);
return
function i3 = I3(x)
i3=besseli(3,x);
return
function i4 = I4(x)
i4=besseli(4,x);
return

% % (unmodified) bessel functions of first kind
% function i0 = I0(x)
% i0=besselj(0,x);
% return
% function i1 = I1(x)
% i1=besselj(1,x);
% return
% function i2 = I2(x)
% i2=besselj(2,x);
% return
% function i3 = I3(x)
% i3=besselj(3,x);
% return
% function i4 = I4(x)
% i4=besselj(4,x);
% return






