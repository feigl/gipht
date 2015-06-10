function j = loopfun(i,n)
%function j = loopfun(i,n)
% like mod(i,n) but return n instead of 0
j=mod(i,n);
if j==0
   j=n;
end
return

