function lr=logproposal2(m1,m2,mstep)
%function lr=logproposal2(m1,m2,mstep)
% 
% given two sets of model parameters m1 and m2, return the log of their
% relative likelihood
% 
% from Example 11.4
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% lr=logr(x,y)
%
% For this problem, we'll use a multivariate normal generator, with 
% standard deviations specified by the vector step.

% Note that logproposal.m and generate.m
% are closely tied to each other.
% 
% 20140107 modified by Kurt Feig
%
lr=(-1/2)*sum((m1-m2).^2./mstep.^2);

return
end
