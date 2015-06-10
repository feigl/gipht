function nerr = check_struct(S,verbose)
% verify that S is an acceptable structure
if exist('verbose','var') ~= 1
    verbose = 0;
end
nerr = 0;
if isstruct(S) ~= 1
    warning('Not a structure');
    nerr = nerr+1;
else
    fn=fieldnames(S);
    if verbose == 1; fprintf(1,'Field    nr nc n(NaN) min max mean std\n'); end;

    for i=1:numel(fn)
        F1=getfield(S,fn{i});
        if isnumeric(F1) == 1
            if numel(find(isfinite(F1)==0)) > 0
                % these are not defined
                if numel(strfind('p1',char(fn{i})))==0 ...
                    && numel(strfind('sigma',char(fn{i})))==0 ...
                    && numel(strfind('mx',char(fn{i})))==0 ...
                    && numel(strfind('my',char(fn{i})))==0 ...
                    && numel(strfind('mz',char(fn{i})))==0                   
                    warning(sprintf('Found NaN in field %s\n',char(fn{i})));
                end
                nerr = nerr + 1;
            end
            if verbose == 1
                fprintf(1,'%3d %20s %5d %5d %8d %#14.6e %#14.6e %#14.6e %#14.6e\n'...
                    ,i,char(fn{i}),size(F1,1),size(F1,2),numel(find(isfinite(F1)==0))...
                    ,nanmin(F1),nanmax(F1),nanmean(F1),nanstd(F1));
            end
        end
    end
end
return
