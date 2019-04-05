function [ F_calc, df1, df2, F_alpha] = Ftest_eqvar( var1, var2, n, alpha)
%[ F_calc, df1, df2, F_alpha] = Ftest_eqvar( var1, var2, n, alpha)
%   Test whehter model 2 provides a significantly better fit than model 1
%   based on variance of residuals
%Elena Baluyut, 11-17-14
% edited ECR 20180613

    disp('Test for equality of variances:')
    disp('Ho = two normal populations have the same variance')
    disp('Ha = second normal population has higher variance')
    
    % define degrees of freedom
    df1 = n-1;
    df2 = n-1;

    % compute F test statistics
    F_calc = var2/var1; 

    F_alpha = icdf('f', 1-alpha, df1, df2);
    
    if F_calc > F_alpha
        disp('reject Ho, model 2 has significantly higher variance than model 1')
    else
        disp('cannot reject Ho at this alpha, model 2 does not have significantly higher variance than model 1 at this alpha')
    end

% end

return
