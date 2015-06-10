function [U2, teststring] = vonmisesness2(theta1)
%function [U2, teststring] = vonmisesness2(theta)
% Test if directional theta theta in radians follow Von Mises distribution
%
%
% The null hypothesis is that theta is distributed as Von Mises
%
% Page 84 of
% N. I. Fisher, Statistical Analysis of Circular Data, 
% Cambridge U. Press (1993)
%
%
% 
% Example 1
%
% x = 0.2*pi*randn([1000 1]);
% th = angle(complex(cos(x),sin(x)));
% [U2,String]  = vonmisesness2(th)
% 
% Example 2
%
% th = randraw('vonmises', [pi/2, 3], [1 1e5]);
% [U2,String]  = vonmisesness2(th)
% 
% 
% Kurt Feigl 2008-MAR-30
% 2009-JAN-12 prune NaN

format long
U2 = NaN;
teststring = 'undefined';

% prune NaN
theta = reshape(theta1,numel(theta1),1);
iok = find(isfinite(theta1)==1);
theta = theta1(iok);

% % prune
% if isreal(theta1) == 0
%     warning('data must be real');
%     return
% end
% 2011-NOV-30 take real part
if isreal(theta) ~= 1
    warning(sprintf('Found complex values in theta. Max abs(imaginary) is %g\n',max(imag(theta))));
    theta = real(theta);
end


% % example of ants
% theta = zeros(10,10);
% theta( 1,:) = [330 290  60 200 200 180 280 220 190 180]';
% theta( 2,:) = [180 160 280 180 170 190 180 140 150 150]';
% theta( 3,:) = [160 200 190 250 180  30 200 180 200 350]';
% theta( 4,:) = [200 180 120 200 210 130  30 210 200 230]';
% theta( 5,:) = [180 160 210 190 180 230  50 150 210 180]';
% theta( 6,:) = [190 210 220 200  60 260 110 180 220 170]';
% theta( 7,:) = [ 10 220 180 210 170  90 160 180 170 200]';
% theta( 8,:) = [160 180 120 150 300 190 220 160  70 190]';
% theta( 9,:) = [110 270 180 200 180 140 360 150 160 170]';
% theta(10,:) = [140  40 300  80 210 200 170 200 210 190]';
% convert to column vector
%theta = reshape(theta,numel(theta),1);
% convert to radians from degrees
%theta = pi * theta / 180;

% sample size
n = numel(theta);
if n < 50
    teststring = sprintf('sample size n (%d) is too small. Need n greater than 50\n',n);
    warning(teststring);
    return
end

% mean direction in radians
md = mean_direction(theta);
%fprintf(1,'mean direction in degrees %f\n',md*180/pi)

% mean resultant length
Rbar = mean_resultant_length(theta);  

% concentration parameter in Von Mises Distribution
kappa = batschelet_inv(Rbar);

%kappa = 1.54;  % adjusted value for example

% shift mean to origin
theta2 = angle(complex(cos(theta - md),sin(theta - md)));

%  Frequencies (4.33)
F=zeros(n,1);
for i=1:n
  F(i) = von_mises_cdf(theta2(i),  0, kappa);
  %fprintf(1,'%12.4e %12.4e\n',theta2(i),F(i));
end

% sort frequencies into increasing order (4.34)
z = sort(F);
%figure;plot(z,'+k');

% calculate test statistic (4.35)
sum = 0;
zbar = mean(z);
for i=1:n
   sum = sum + ( z(i) - (2*i - 1)/(2*n) )^2;
end
U2 = sum - n*(zbar-0.5)^2 + 1/(12*n);

%critical values at significance level alpha = 0.05
%Case I: mu unknown, kappa known
% if kappa < 0.25
%    crit = 0.133;
% elseif kappa < 0.5
%    crit = 0.135;
% elseif kappa < 1.0
%    crit = 0.139;
% elseif kappa < 1.5
%    crit = 0.144;
% elseif kappa < 2.0
%    crit = 0.147;  
% elseif kappa < 4.0
%    crit = 0.153;
% elseif kappa < 1.0E99
%    crit = 0.157;
% else
%    error(sprintf('Unknown kappa %e\n',kappa));
% end
%Case III: mu unknown, kappa unknown
if kappa < 0.25
   crit = 0.061;
elseif kappa < 0.5
   crit = 0.066;
elseif kappa < 1.0
   crit = 0.079;
elseif kappa < 1.5
   crit = 0.092;
elseif kappa < 2.0
   crit = 0.101;  
elseif kappa < 4.0
   crit = 0.113;
elseif kappa < 1.0E99
   crit = 0.117;
else
   warning(sprintf('Unknown kappa %e Setting crit to NaN\n',kappa));
   crit = NaN;
end
   

if U2 < crit
   teststring = sprintf('Data ARE compatible with Von Mises distribution because U2 = %g is less than critical value %10.4f\n',U2,crit);
else
   teststring = sprintf('Data are NOT compatible with Von Mises distribution because U2 = %g is greater than critical value %10.4f\n',U2,crit);
end

if nargout < 2
   fprintf(1,'%s\n',teststring);
end

return






