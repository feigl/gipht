function [ F_calc, dfn, dfd, F_alpha] = Ftest( var1, var2, n1, n2, alpha)
%Ftest_regression:  F-test for regression analysis
%   Test whether the fit of model 2 is better than the fit of model 1,
%   where model 2 is larger.  SSWR is sum of squared weighted residuals, p
%   is number of parameters in each model (tbreaks), and n is number
%   of data points, alpha is the level of test such that F_calc falls 
%in the 1-alpha region of the F distribution (e.g., .100, .050, .025, .005)
%Elena Baluyut, 11-17-14
disp('Performing regression F test on models 1 and 2, Assume Model 1 is nested in Model 2.')
disp('Ho = no significant difference in fits')
disp('Ha = Model 2 provides a better fit than Model 1')
F_calc = var1/var2;
dfn = n1 - 1;
dfd = n2 - 1;
alpha = alpha;
F_alpha = icdf('f', 1-alpha, dfn, dfd);
if F_calc > F_alpha
    disp('reject Ho, model 2 adds significant fit')
else
    disp('cannot reject Ho at this alpha, not a significant amount of fit contributed by model 2 at this alpha')
end

