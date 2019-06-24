function [ F_calc, df1, df2, F_alpha, H] = Ftest_modelcomp( SSWR1, SSWR2, p1, p2, n, alpha)
%   [F_calc, df1, df2, F_alpha, H] = Ftest_modelcomp( SSWR1, SSWR2, p1, p2, n, alpha)
%   Test whether the fit of model 2 is better than the fit of model 1,
%   where model 2 is larger (more complex).  SSWR is sum of squared weighted residuals, p
%   is number of parameters in each model, and n is number
%   of data points, alpha is the level of test such that F_calc falls 
%in the 1-alpha region of the F distribution (e.g., .100, .050, .025, .005)
%Elena Baluyut, 11-17-14
% edited ECR 20180613


    disp('Performing regression F test on models 1 and 2, Assume Model 1 is less complex than Model 2.')
    disp('Ho = Model 2 does not provide significantly better fit than Model 1')
    disp('Ha = Model 2 provides a better fit than Model 1')

    % define degrees of freedom
    df1 = n-p1;
    df2 = n-p2;
    
    % compute F test statistics
    F_calc = ((SSWR1 - SSWR2)/(df1-df2))/((SSWR2)/(df2));
    
    F_alpha = icdf('f', 1-alpha, df1, df2);

    if F_calc > F_alpha
        disp('reject Ho, model 2 adds significant fit')
        H = 1;
    else
        disp('cannot reject Ho at this alpha, not a significant amount of fit contributed by model 2 at this alpha')
        H = 0;
    end

return
