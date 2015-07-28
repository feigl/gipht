function [ F_calc, dfn, dfd, F_alpha] = Ftest_nl( SSWR1, SSWR2, p1, p2, n, alpha)
%Ftest_regression:  F-test for regression analysis
%   Test whether the fit of model 2 is better than the fit of model 1,
%   where model 2 is larger (more complex).  SSWR is sum of squared weighted residuals, p
%   is number of parameters in each model, and n is number
%   of data points, alpha is the level of test such that F_calc falls 
%in the 1-alpha region of the F distribution (e.g., .100, .050, .025, .005)
%Elena Baluyut, 11-17-14
disp('Performing regression F test on models 1 and 2, Assume Model 1 is less complex than Model 2.')
disp('Ho = no significant difference in fits')
disp('Ha = Model 2 provides a better fit than Model 1')

dfn = n-p1;
dfd = n-p2;

if p1 == p2
    F_calc = SSWR1/SSWR2; 
else
    F_calc = ((SSWR1 - SSWR2)/(p2 - p1))/((SSWR2)/(n-p2));
end

alpha = alpha;

F_alpha = icdf('f', 1-alpha, dfn, dfd);

if F_calc > F_alpha
    disp('reject Ho, model 2 adds significant fit')
else
    disp('cannot reject Ho at this alpha, not a significant amount of fit contributed by model 2 at this alpha')
end


