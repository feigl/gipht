function [val,irow,jcol] = find_matrix_value_indices(A,fn)
% return maximum value of 2-D array A and its indices
%https://www.mathworks.com/matlabcentral/answers/63247-get-max-value-and-index-of-multidimensional-array
if numel(size(A)) == 2
    [n,m] = size(A);
    switch fn
        case {'min','nanmin','max','nanmax'}         
            [val,jcol] = feval(fn,feval(fn,A));
            [~,irow] = feval(fn,A(:,jcol));
        case {'mean','nanmean','std','nanstd','median','nanmedian','mode','nanmode'}         
            val = feval(fn,reshape(A,n*m,1));
            irow=nan;
            jcol=nan;
        otherwise
            error('unknown function %s\n',fn);
    end
else
    error('expecting 2-D array');
end
return
end

