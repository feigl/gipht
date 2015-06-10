function DST1 = bootstrap_dst(DST0)
% Given a DST structure with n elements in each field, 
% return another DDST0T structure with a bootstrap resampling
% with replacement
% the number of samples is the same
% 2012-JAN-03 Kurt Feigl
verbose = 0;
DST1 = struct('dummy',NaN);
if isstruct(DST0) ~= 1
    warning('Not a structure');
    nerr = nerr+1;
else
    %DST1 = orderfields(DST0);
    fn=fieldnames(DST0);
    if verbose == 1
        fprintf(1,'Field    nr nc n(NaN) min max mean std\n');
    end

    % define pointers to elements
    F1=getfield(DST0,fn{1});
    ndat = numel(F1);
    jj = ceil(ndat*rand(size(F1)));

    % apply same pointers to all fields
    for i=1:numel(fn)
        F1=getfield(DST0,fn{i});
        clear F2;
        if isnumeric(F1) == 1
            if numel(find(isfinite(F1)==0)) > 0
                % these are not defined
                if numel(strfind('p1',char(fn{i})))==0 && numel(strfind('sigma',char(fn{i})))==0
                    warning(sprintf('Found NaN in field %s\n',char(fn{i})));
                end
                nerr = nerr + 1;
            end
            if verbose == 1
                fprintf(1,'%3d %20s %5d %5d %8d %#14.6e %#14.6e %#14.6e %#14.6e\n'...
                    ,i,char(fn{i}),size(F1,1),size(F1,2),numel(find(isfinite(F1)==0))...
                    ,nanmin(F1),nanmax(F1),nanmean(F1),nanstd(F1));
            end
            if numel(F1) == ndat
                % choose new values
                F2 = colvec(F1(jj));
                % set the field in the new structure
                DST1 = setfield(DST1,fn{i},F2); 
                if verbose == 1
                    fprintf(1,'Getting field named %s\n',fn{i});
                end
                F3=getfield(DST1,fn{i});
                if verbose == 1
                    fprintf(1,'%3d %20s %5d %5d %8d %#14.6e %#14.6e %#14.6e %#14.6e\n'...
                        ,i,char(fn{i}),size(F3,1),size(F3,2),numel(find(isfinite(F3)==0))...
                        ,nanmin(F3),nanmax(F3),nanmean(F3),nanstd(F3));
                end
                else
                    warning(sprintf('Mismatch in number of fields. Found %d expected %d\n',numel(F1),ndat));
                end
       end
    end
end
return
