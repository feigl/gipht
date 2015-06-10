function ix=regexpindex(strs,pat)
%function ix=regexpindex(strs,pat)
% return indices of cells in strs matching regexp pat
% Kurt Feigl 2009-JUL-27
h=regexp(strs,pat);
k=0;
for i=1:numel(h)
   if numel(h{i}) > 0
      k=k+1;
      ix(k) = i;
   end
end
if k==0
   ix=[];
end
return
