function normalized_residual=fun_normalized_residuals_with_recording(params, PST, DST, TST)
% return residuals
% for use with LSQNONLIN
% 20200609 Kurt Feigl
% keep values in memory
%persistent funit
%persistent fnameout
persistent ncallcount

% persistent param_history 
% persistent cost1 
% persistent nchi1

fnameout = PST.fnameout;

% count the number of calls to this function
if isempty(ncallcount) == true
    ncallcount = 0;
else
    ncallcount = ncallcount + 1;
end

% count parameters
if numel(params) ~= numel(PST.p1)
    error('miscount of parameters');
else
    mparams = numel(PST.p1);
end

if ncallcount == 0   
    for ip=1:mparams
        fprintf(1,' %-8s',PST.name{ip});
    end
    fprintf(1,'\n');
else
    for ip=1:mparams
        %fprintf(1,'    Re-setting %d %-8s from %12.4e to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),params(ip),PST.description{ip});
        %fprintf(1,'    Re-setting %d %-8s from %12.4e to %12.4e\n',ip,PST.name{ip},PST.p1(ip),params(ip));
        fprintf(1,' %12.4e',params(ip));
    end
        fprintf(1,'\n');
end

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
modeled = Hfitfun(PST,DST,TST);
modeled = reshape(modeled,ndata,1);

% calculate residual
residual = DST.obs-modeled;

% normalize residual by measurement uncertainty, neglecting correlation
normalized_residual = residual ./ DST.sig;

% sum of squares of residuals
sumsqrnres = normalized_residual' * normalized_residual;

% misfit
sqrtnchi2 = sqrt(sumsqrnres/ndata);


%% open file 
if ncallcount == 0 
    % open file, flush, and write header
    funit = fopen(fnameout,'wt');
    if funit ~= -1
        fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        fprintf(funit,'NCALLCOUNT,SUMSQRNRES,SQRTNCHI2');
        for i=1:numel(PST.name)
            fprintf(funit,',%s',PST.name{i});
        end
        fprintf(funit,'\n');
        fclose(funit);
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
else
    % open file to append
    funit = fopen(fnameout,'At');
    if funit ~= -1
        %fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        if numel(isfinite(PST.p1)) == mparams && isfinite(sumsqrnres) == true && isfinite(sqrtnchi2) == true
            fprintf(funit,'%09d',ncallcount);
            fprintf(funit,',%12.7E',sumsqrnres);
            fprintf(funit,',%12.7E',sqrtnchi2);
            for i=1:numel(PST.p1)
                fprintf(funit,',%12.7E',PST.p1(i));
            end
            fprintf(funit,'\n');
        end
        fclose(funit);    
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
end
return
end

