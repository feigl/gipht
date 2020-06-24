
ndata1 = 10;
ndata2 = 10;
ndof1 = 9;
ndof2 = 9;

var1 = 1.0;
% var2 = 2.0;
var2 = 1.0;

variance_ratio = var2/var1

X1 = (var1^2)*randn(ndata1,1);
X2 = (var2^2)*randn(ndata2,1);

alpha = 1-0.682

% alternative hypothesis
% 'right' -- "variance of X is greater than variance of Y"
%                        (right-tailed test)
[H,P,CI,STATS] = vartest2(X2,X1,'alpha',alpha,'tail','right')
H

fcrit1 = ftest_critical(alpha,variance_ratio,ndof1,ndof2)
fcrit2 = icdf('F',1-alpha,ndof1,ndof2)



pvalT1 = fpval(variance_ratio, ndof1, ndof2)
pvalT2 = 1. - cdf('F',variance_ratio, ndof1, ndof2)
pvalCrit = cdf('F',fcrit2,ndof1,ndof2)





