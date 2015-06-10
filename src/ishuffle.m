function j=ishuffle(i)
% given an array of non-zero positive integers i,
% shuffle their order randomly and return in j
r=randn(size(i));
[r,k]=sort(r,1,'ascend');
j=i(k);
return

