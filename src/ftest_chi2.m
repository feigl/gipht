function [Fcritical,H,test_string] = ftest_chi2(alpha,chiSquare1,chiSquare2,ndof1,ndof2)
%% given two chi-square statistics, perform one-sided, two-sample F-test of null hypothesis:
%    H0: chiSquare1 equals chiSquare2
% versus alternative hypothesis:
%    H1: chisquare1 < chiSquare2
%
%function [fcritical,H] = ftest_chi2(alpha,chiSquare1,chiSquare2,ndof1,ndof2)
%
% Inputs:
%  alpha      == level of significance (e.g., 0.05 for 95% confidence)
%  chiSquare1 == chi-squared statistic for 1st sample
%  chiSquare2 == chi-squared statistic for 2nd sample 
%    for clarity, assume let chiSquare2 > chiSquare1
%    chiSquare is defined by:
%    Equation (15.1.5) Press et al. (1992)
%    Equation (11.3) Bevington and Robinson (2003)
%    chiSquare = sum(ydobs-ymod)./sig)^2 
%    chiSquare = nres' * nres; where nres = ydobs-ymod)./sig
%  ndof1      == number of degrees of freedom (N-M) for 1st sample
%  ndof2      == number of degrees of freedom (N-M) for 2nd sample
%
% Outputs:
%  fcritical  == critical value of F statistic to reject H0
%  H          == 0 if null hypothesis H0 fails to be rejected,
%             == 1 if null hypothesis H0 is rejected in favor of alternative hypothesis H1
%
% References:
%
% Bevington, P. R., and D. K. Robinson (2003), Data reduction and error analysis for the physical sciences, 3rd ed., xi,
% 320 p. pp., McGraw-Hill, Boston.  http://www.loc.gov/catdir/description/mh024/2002070896.html
%
% Press, W. H., S. A. Teukolsky, B. P. Flannery, and W. T. Vetterling (1992), Numerical recipes in Fortran 77: volume 1,
% volume 1 of Fortran numerical recipes: the art of scientific computing, Cambridge university press.
%
% Larsen, R. J., and M. L. Marx (1986), An introduction to mathematical statistics and its applications, 2nd ed., x, 630
% p. pp., Prentice-Hall, Englewood Cliffs, N.J.
%
% 
% 20200702 Kurt Feigl

%% initialize
narginchk(5,5);
Fcritical = nan;
H=nan;
test_string = '';

%% check inputs
if alpha >= 1.0 || alpha <= 0.
    warning('signficance level alpha must be strictly between 0. and 1.0');
    return
else
    confidence = 100. * (1.0 - alpha);
end
if ndof1 < 1 || ndof2 < 1
    warning('number of degrees of freedom must be at least 1');
    return
end
if chiSquare1 <= 0 || chiSquare2 <= 0
    warning('chiSquare statistic must be greater than 0');
    return
end
if chiSquare1 > chiSquare2 
    warning('One-sided test requires chiSquare1 <= ChiSquare2');
    return
end

%% calculate observed values of f
% Bevington and Robinson (2003) write: If two statistic chiSquare1 and chiSquare2 which follow the chi-square
% distribution have been determined, the ratio of the reduced chi-squareds, rchiSquare1 and rchiSquare2, is distributed
% according to the F distribution:
% f = rchiSquare1/rchiSquare2
F12 = (chiSquare1/ndof1) / (chiSquare2/ndof2);
F21 = (chiSquare2/ndof2) / (chiSquare1/ndof1);

%% find probability of observed values occuring
PF12 = cdf('F',F12, ndof1, ndof2);
PF21 = cdf('F',F21, ndof1, ndof2);

% critical value
Fcritical = icdf('F',1.0-alpha, ndof1, ndof2);

% verify
PFcritical = cdf('F',Fcritical,ndof1,ndof2);

discrepancy = PFcritical - (1.0 - alpha);
if abs(discrepancy) > 10*eps
    warning(sprintf('Discrepancy %12.4E is large.\n',discrepancy));
end

fprintf(1,'Performing one-sided, two-sample F-test of null hypothesis:\n');
fprintf(1,'   H0:          chiSquare1 %10.2g equals         chiSquare2 %10.2g\n',chiSquare1,chiSquare2);
fprintf(1,'   H0: redcuced chiSquare1 %10.2g equals reduced chiSquare2 %10.2g\n',chiSquare1/ndof1,chiSquare2/ndof2);
fprintf(1,'versus alternative hypothesis:\n');
fprintf(1,'   H1: chisquare1 < chiSquare2\n');

% two-sided test
% if PF12 <= alpha/2. || PF21 >= (1. - alpha/2)
% one-sided test
if PF21 >= (1. - alpha)
    H = 1;
    test_string = sprintf('Null hypothesis rejected with %.0f %% confidence',confidence);
else
    H = 0;
    test_string = sprintf('Null hypothesis fails to be rejected with %.0f %% confidence',confidence);
end
fprintf(1,'%s\n',test_string);


%% set up limits for plot
%fvec=logspace(log10(nanmin([F12,F21])/2.),log10(nanmax([F12,F21])*2.),100);
Fmin = nanmin([F12,F21,Fcritical])/2.;
Fmax = nanmax([F12,F21,Fcritical])*2.;
fvec=logspace(log10(Fmin),log10(Fmax),100);
Pc=cdf('F',fvec,ndof1,ndof2);
Pd=pdf('F',fvec,ndof1,ndof2);
Pmin=nanmin([Pc,Pd]);



figure;
%hold on;
loglog(fvec,Pc,'k-','LineWidth',2);
hold on;
loglog(fvec,Pd,'b-','LineWidth',2);
loglog([Fcritical, Fcritical],[Pmin, PFcritical],'r--','LineWidth',3);
%loglog([F12, F12],[1.E-6, PF12],'k--','LineWidth',2);
loglog([F21, F21],[Pmin, PF21],'m:', 'LineWidth',3);
xlabel('f');
ylabel('P');
set(gca,'FontWeight','bold','FontSize',12);
title(sprintf('F(%d, %d) %s',ndof1, ndof2,test_string));
%legend('CDF','PDF','critical','F12','F21');
legend('CDF','PDF','critical','F21','Location','SouthWest');
%legend('PDF','critical','F12','F21');

print(gcf,'-dpdf',sprintf('%s.pdf',mfilename),'-r600','-fillpage','-painters');

return
end
