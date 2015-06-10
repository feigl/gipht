function pnames2 = truncate_parameter_names(pnames1)
% function pnames2 = truncate_parameter_names(pnames1)
%   Truncate and pad list of parameter names
mparam = numel(pnames1);

% change spaces to underscores
for i=1:mparam
    tmpstr = pnames1{i};
    nchar = numel(tmpstr);
    if nchar < 32
       tmpstr(nchar+1:32) = ' '; % pad to 32 characters
    else
       if nchar > 32
           warning(sprintf('Truncating %d characters from parameter name %s\n',nchar-32,tmpstr));
       end
       tmpstr = tmpstr(1:32); % truncate to 32 characters
    end
    tmpstr2 = strrep(tmpstr,' ','_');
    pnames2{i} = tmpstr2;
end

return

