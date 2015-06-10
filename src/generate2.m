function mtrial=generate2(mcurrent,msig)
%function mtrial=generate2(mcurrent,msig)
% generate trial values of model parameters, given uncertainties msig
% from
% Example 11.4
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% y=generate(x)
%
%
% For this problem, we'll use a multivariate normal generator, with 
% standard deviations specified by the vector msig.
%
% Note that logproposal.m and generate.m
% are closely tied to each other.
%
% modified 20140107 Kurt Feigl

mtrial=mcurrent+msig.*randn(numel(msig),1);
return
end
