function s=lststr(t)
%function s=lststr(t)
% find the last string s starting with a blank in t
% last updated 2013-05-20 Kurt


m = length(t);
if m > 0
    
    % last item is a white space
    if isspace(t(m))
        while isspace(t(m)),
            m = m -1;
        end
    end;
    
    while ~ isspace(t(m)),
        m = m -1;
    end
    
    s = t(m:length(t));
else
    s = [];
end

return;
