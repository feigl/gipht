function [ T_calc, df, T_alpha, Sp] = Ttest_2tail( m1, m2, Do, S1, S2, n1, n2, alpha)
%Ttest_2tail performs T test for 2 sided distribution
%   m is the slope of the interval
%   Do is the difference tested between slopes 
%   S is the standard deviation or uncertainty
%   n is the number of data points per invertal 
%   alpha is confidence level; (1-alpha)% confidence
% Original: Elena C. Reinisch
% 2020/01/05: Change disp to fprintf

fprintf(1,'Performing regression T test on slopes for current tbreaks with significance level alpha = %.4f\n',alpha);
fprintf(1,'Ho == Slopes are equal.\n')
fprintf(1,'Ha == Difference in slopes is significantly different from zero.\n')
Sp = sqrt(((n1-1)*S1^2+(n2-1)*S2^2)/(n1+n2-2));
T_calc = (m1-m2-Do)/(Sp*sqrt(1/n1+1/n2));
df = n1+n2-2;
T_alpha = icdf('t', 1-alpha/2, df);
if abs(T_calc) > T_alpha
    fprintf(1,'Reject Ho, significant difference in slopes.\n')
else
    fprintf(1,'Cannot reject Ho at this alpha, not a significant difference in slopes at this alpha.\n')
end
return
end