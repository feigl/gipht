function ix=regexpiindex(strs,pat)
%function ix=regexpiindex(strs,pat)
% return indices of cells in strs matching regexpi pat 
% Thi version is insensitive to case
% Kurt Feigl 20151122
h=regexpi(strs,pat);
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
