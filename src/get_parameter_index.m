function pindex = get_parameter_index(pname1,pnames)
% return the index of parameter pname1 in list pnames

% now look for index or index array
%pindex = strmatch(pname1,pnames);
%pindex = regexpi(pnames,pname1,'start');
%pindex=regexpindex(pnames,'\w*epoch')

if numel(pname1) > 32
    warning(sprintf('Number of characters %d in %s\n',numel(pname1),pname1));
    pname1 = pname1(1:32);
end
        
pindex=regexpindex(pnames,pname1);
% if numel(pindex) == 0
%     warning(sprintf('Parameter name not found: %s\n',pname1));
% end

return






