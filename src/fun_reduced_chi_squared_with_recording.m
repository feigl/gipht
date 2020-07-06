function rchi2 = fun_reduced_chi_squared_with_recording(params, PST, DST, TST)
% return reduced chi squared statistic
% for use with FMINCON
%
% inputs:
%    params   == current values of model parameters
%    PST      == structure containing model parameters
%    DST      == structure containing data and metadata
%    TST      == structure containing temporary items
%
% chiSquare is defined by:
%    Equation (15.1.5) Press et al. (1992)
%    Equation (11.3) Bevington and Robinson (2003)
%    chiSquare = sum(ydobs-ymod)./sig)^2 
%    chiSquare = nres' * nres; where nres = ydobs-ymod)./sig
% reduced chi-squared is chiSquare/ndof
%    
% References:
%
% Bevington, P. R., and D. K. Robinson (2003), Data reduction and error analysis for the physical sciences, 3rd ed., xi,
% 320 p. pp., McGraw-Hill, Boston.  http://www.loc.gov/catdir/description/mh024/2002070896.html
%
% Press, W. H., S. A. Teukolsky, B. P. Flannery, and W. T. Vetterling (1992), Numerical recipes in Fortran 77: volume 1,
% volume 1 of Fortran numerical recipes: the art of scientific computing, Cambridge university press.
% 
% 20200705 Kurt Feigl
% keep values in memory
persistent ncallcount
persistent ticTotal;

fnameout = PST.fnameout;

% remember when we started
if isempty(ticTotal) == true
    ticTotal = tic;
else
    tElapsed = toc(ticTotal);
    tElapsedFormatted = seconds(tElapsed);
    tElapsedFormatted.Format = 'dd:hh:mm:ss.SSS';
end

% count the number of calls to this function
if isempty(ncallcount) == true
    ncallcount = 0;
else
    ncallcount = ncallcount + 1;
end

% count parameters
if numel(params) ~= numel(PST.name)
    ncallcount
    params
    numel(params)
    PST.name
    error('miscount of parameters');
else
    mparams = numel(PST.p1);
    
end

% if ncallcount == 0
%     for ip=1:mparams
%         fprintf(1,' %-8s',PST.name{ip});
%     end
%     fprintf(1,'\n');
% else
%     for ip=1:mparams
%         %fprintf(1,'    Re-setting %d %-8s from %12.4e to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),params(ip),PST.description{ip});
%         %fprintf(1,'    Re-setting %d %-8s from %12.4e to %12.4e\n',ip,PST.name{ip},PST.p1(ip),params(ip));
%         fprintf(1,' %12.4e',params(ip));
%     end
%         fprintf(1,'\n');
% end

% copy current value of parameters and rescale to values with units
PST.p1 =  PST.p0 + params .* PST.scale;

% get a handle on fitting function
Hfitfun = PST.Hfitfun;
if isa(Hfitfun,'function_handle') == false
    fprintf(1,'Making function handle\n');
    Hfitfun = @PST.Hfitfun;
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

% sum of squares of normalized residuals
chi2 = normalized_residual' * normalized_residual;

% number of degrees of freedom
ndof = ndata - mparams;

% reduced chi squared
rchi2 = chi2/ndof;

%% open file 
if ncallcount == 0 
    % open file, flush, and write header
    funit = fopen(fnameout,'wt');
    if funit ~= -1
        fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        fprintf(funit,'ElapsedSeconds, nCallCount, rchi2');
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
        funits = [1, funit];
        %fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        if numel(isfinite(PST.p1)) == mparams && isfinite(rchi2) == true
            for j=1:numel(funits)
                funit1 = funits(j);
                %fprintf(funit,'%s',datestr(now,30));
                %fprintf(funit,',%s',tElapsedFormatted);
                fprintf(funit1,'%12.3E',tElapsed);
                fprintf(funit1,',%09d',ncallcount);
                fprintf(funit1,',%12.7E',chi2);
                for i=1:numel(PST.p1)
                    fprintf(funit1,',%12.7E',PST.p1(i));
                end
                fprintf(funit1,'\n');
            end
        end
        fclose(funit);    
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
end
return
end

