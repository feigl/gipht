function z = zscore(x, s)
%function z = zscore(x, s)
% given random variable x with uncertainty s, return Z-score
% 20200402 Kurt Feigl
% example 1:
% x = 100*randn(100,1)+50;
% z = zscore(x);
% figure
% subplot(2,1,1);
% histogram(x);
% subplot(2,1,2);
% histogram(z);
% example 2:
% X = randn(10)
% S = nanstd(colvec(X))*ones(size(X))
% Z=zscore(X,S)
if nargin == 1
    s = nanstd(colvec(x));
end
z = (x - nanmean(colvec(x))) ./ s;
return
end

