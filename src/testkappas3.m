 function P = testkappas3(Rbar1, n1, Rbar2, n2)
%function P = testkappas3(Rbar1, n1, Rbar2, n2)
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
%               abs(eta) > 1.96
% 
% This routine handles:
%  Case I   for 0    <  Rbar <  0.45 i.e. 0 < kappa < 1
%  Case II  for 0.45 <= Rbar <= 0.70 i.e. 1 < kappa < ?
%  Case III for 0.7  >  Rbar   
%
% Kurt Feigl 2010-JAN-11
%
% Test Case for n1 = n2 = 10
%    P=testkappas3(0.2,10,0.1,10)
%    eta = icdf('normal',P,0,1)
%
%  n=10;
%  Rbar1=0.2;
%  Rbar2=0.1;
%  eta1 = sqrt(2*(n-4)/3) * (asin(Rbar1*sqrt(3/2))-asin(Rbar2*sqrt(3/2)));
%  eta1 =
% 
%     0.2494
%
% so eta = eta1 QED

% Test assumptions
if nargin ~= 4
    warning('Need 4 arguments');
    P = NaN;
    return;
end
if n1 < 5 || n2 < 5
    warning('Sample size too small.');
    P = NaN;
    return;
end
if Rbar1 > 1.0 ||  Rbar1 < 0.0 
    P = NaN;
    warning(sprintf('Undefined Rbar1 %12.4e\n',Rbar1));
    return
end
if Rbar2 > 1.0 ||  Rbar2 < 0.0 
    P = NaN;
    warning(sprintf('Undefined Rbar2 %12.4e\n',Rbar2));
    return
end



% Case III
if  Rbar1 > 0.7 || Rbar2 > 0.70
       % numerator
       top = (n1 - Rbar1)/(n1-1);
       
       % denominator
       bot = (n2 - Rbar2)/(n2-1);
       
       % ratio 
       fstat = top/bot;
       P = cdf('F',fstat,n1-1,n2-1);
% case I
elseif (Rbar1 < 0.45 && Rbar2 < 0.45)  || mean([Rbar1 Rbar2]) < 0.45      
   % numerator
   top = 2 * ( g1(2*Rbar1) - g1(2*Rbar2));

   % denominator
   bot = sqrt(3 * (1/(n1 - 4) + 1/(n2 - 4)));

   % ratio
   eta = top/bot;
   
   % probability
   if isreal(eta) == 1
       P = cdf('normal',eta,0,1);
   else 
       P = NaN;
   end
%elseif Rbar1 >= 0.45 & Rbar2 >= 0.45 & Rbar1 <= 0.70 & Rbar2 <= 0.70 | mean([Rbar1 Rbar2]) <= 0.70
else
   % numerator
   top =  g2(Rbar1) - g2(Rbar2);

   % denominator
   c3 = 0.893;
   bot = c3*sqrt(1/(n1 - 3) + 1/(n2 - 3));

   % ratio
   eta = top/bot;
%   eta = NaN;
   P = cdf('normal',eta,0,1);
end


%fprintf(1,'Rbar1 = %12.4e P = %12.4e\n',Rbar1,P);
return



function varstab1= g1(R)
% Variance Stabilizing Transformation
% Mardia and Jupp, Directional Statistics, Wiley, 1999
% equation 4.8.40 p. 81
% Kurt Feigl 2008-JAN-10
a = sqrt(3/8);
varstab1 = asin(a*R);
return

function varstab2= g2(R)
% Variance Stabilizing Transformation
% Mardia and Jupp, Directional Statistics, Wiley, 1999
% equation 4.8.40 p. 82
% Kurt Feigl 2008-FEB-11
c1 = 1.089;
c2 = 0.258;
t = (R-c1)/c2;
varstab2 = asinh(t);
return




