function normalized_residual=fun_normalized_residuals(params, PST, DST, TST)
% return residuals
% for use with LSQNONLIN
% 20200609 Kurt Feigl
% keep values in memory
%persistent funit
persistent fnameout
% persistent param_history 
% persistent cost1 
% persistent nchi1

% copy current value of parameters
PST.p1=params;

% get a handle on fitting function
Hfitfun = PST.fitfun;
if isa(Hfitfun,'function_handle') == false
    fprintf(1,'Making function handle\n');
    Hfitfun = @PST.fitfun;
end

% number of data
ndata = numel(DST.obs);

% evaluate fitting function with current values of parameters, make column vector
modeled = reshape(Hfitfun(PST,DST,TST),ndata,1);

% calculate residual
residual = DST.obs-modeled;

% normalize residual by measurement uncertainty, neglecting correlation
normalized_residual = residual ./ DST.sigma;

% sum of squares of residuals
sumsqrnres = normalized_residual' * normalized_residual;

% misfit
sqrtnchi2 = sqrt(sumsqrnres/ndata);

if exist('fnameout','variable') ~= 1
    fnameout = sprintf('%s_out_%s.csv',mfilename,datestr(now,30));
    % open file and flush
    funit = fopen(fnameout,'wt');
    if funit ~= -1
        fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        fprintf(funit,'SUMSQRNRES,SQRTNCHI2,');
        for i=1:numel(PST.name)
            fprintf(funit,'%s,',PST.name{i});
        end
        fprintf(funit,'\n');
        fclose(funit);
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
else
    % open file to append
    funit = fopen(fnameout,'A');
    if funit ~= -1
        %fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        fprintf(funit,'%12.7E,',sumsqrnres);
        fprintf(funit,'%12.7E,',sqrtnchi2);
        for i=1:numel(PST.name)
            fprintf(funit,'%s,',PST.p1{i});
        end
        fprintf(funit,'\n');
        fclose(funit);
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
end
return
end

