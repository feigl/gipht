function normalized_residual=fun_normalized_residuals_comsol(params, PST, DST, TST)
% return residuals
% for use with LSQNONLIN
% 20200609 Kurt Feigl
% keep values in memory
%persistent funit
%persistent fnameout
persistent ncallcount
persistent sumsqrnres
persistent sqrtnchi2
import com.comsol.model.util.*

% persistent param_history 
% persistent cost1 
% persistent nchi1
verbose = 0;

fnameout = PST.fnameout;

% count the number of calls to this function
if isempty(ncallcount) == true
    ncallcount = 0;
else
    ncallcount = ncallcount + 1;
end


% count parameters
if numel(params) ~= numel(PST.name)
    error('miscount of parameters');
else
    mparams = numel(params);
end

if verbose >= 1
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
end

%% load the mph file 
model = mphload(PST.fileNameMPH);

%https://www.comsol.com/forum/thread/128742/avoid-mphload-every-iteration-in-matlab-livelink
% Inside you objective function you can now do of one these:
% 1)
% import com.comsol.model.util.*
% model = ModelUtil.model('MyModel');
% gives you access to the MyModel model on the server using the variable name model in Matlab

%% set parameters
nreset = 0;
for ip=1:mparams  % loop over indices to comsol parameters
    param1 = params(ip);
    if isfinite(param1) == 1
        db = PST.ub(ip)-PST.lb(ip); % range of bounds in gipht parameter
        scale1 = PST.scale(ip);
        if  db > 0
            %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
            pname1 = sprintf('%s',PST.name{ip});
            unit1 = char(model.param.evaluateUnit(pname1));
            desc1 = char(model.param.descr(pname1));
            if verbose >= 2
                %fprintf(1,'    Setting %d %-8s to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),PST.description{ip});
                fprintf(1,'    Setting %d %-8s to %12.4e [%s] %s\n',ip,pname1,params(ip),unit1,desc1);
            end
            %model.param.set(pname1,sprintf('%e [%s]',param1,unit1),desc1);
            model.param.set(pname1,sprintf('%e [%s]',param1*scale1,unit1),desc1);
            nreset = nreset + 1;
        end        
    end
end

%% run comsol solution and save it
%fprintf(1,'%s: starting COMSOL solution with mph file named %s\n',mfilename,PST.fileNameMPH);
%fprintf(1,'%s: starting COMSOL solution...\n',mfilename);
trun0 = tic;
model.sol(PST.soltag).runAll;
mphsave(model);
%fprintf(1,'%s: finished COMSOL solution in %10.1f seconds\n',mfilename,toc(trun0));


%% calculate modeled values at locations of data points
expr = mphgetexpressions(model.param);
sites = DST.sites;
nsites = numel(sites);
ndata = 3 * nsites;
modeledDisplacements = nan(ndata,1);
idata=0;
for isite = 1:nsites
    ix=find(contains(expr(:,1),lower(sprintf('x%4s',sites{isite}))));
    iy=find(contains(expr(:,1),lower(sprintf('y%4s',sites{isite}))));   
    U = mphinterp(model,{'u2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',PST.dset); %east
    V = mphinterp(model,{'v2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',PST.dset); % north component
    W = mphinterp(model,{'w2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',PST.dset); % vertical component
    idata=idata+1;modeledDisplacements(idata) = U;
    idata=idata+1;modeledDisplacements(idata) = V;
    idata=idata+1;modeledDisplacements(idata) = W;
end

if idata ~= ndata
    error('miscount')
end
modeledDisplacements = reshape(modeledDisplacements,ndata,1);

if verbose >= 2
fprintf(1,'Min, Max, Median of displacement in meters %10.3f %10.3f %10.3f\n'...
    ,nanmin(modeledDisplacements) ...
    ,nanmax(modeledDisplacements) ...
    ,nanmedian(modeledDisplacements));
end

if verbose >= 2
    figure;
    histogram(modeledDisplacements);
    xlabel('modeled displacement [m]');
    ylabel('count');
end

% calculate residual
residual = DST.obs -  modeledDisplacements;

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
        fprintf(funit,'NCALLCOUNT, SUMSQRNRES, SQRTNCHI2');
        fprintf(1,    'NCALLCOUNT, SUMSQRNRES, SQRTNCHI2');
        for i=1:numel(PST.name)
            fprintf(funit,',%s',PST.name{i});
            fprintf(1    ,',%s',PST.name{i});
        end
        fprintf(funit,'\n');
        fprintf(1    ,'\n');
        fclose(funit);
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
else
    % open file to append
    funit = fopen(fnameout,'At');
    if funit ~= -1
        %fprintf(1,'Opened file named %s to record parameters.\n',fnameout);
        if numel(isfinite(params)) == mparams && isfinite(sumsqrnres) == true && isfinite(sqrtnchi2) == true
            fprintf(funit,'%09d',ncallcount);
            fprintf(1    ,'%09d',ncallcount);
            fprintf(funit,',%12.7E',sumsqrnres);
            fprintf(1    ,',%12.7E',sumsqrnres);
            fprintf(funit,',%12.7E',sqrtnchi2);          
            fprintf(1    ,',%12.7E',sqrtnchi2);
            for i=1:numel(params)
                % fprintf(funit,',%12.7E',params(i));
                % fprintf(1    ,',%12.7E',params(i));
                % print unscaled values of parameters
                fprintf(funit,',%12.7E',params(i)*PST.scale(i));
                fprintf(1    ,',%12.7E',params(i)*PST.scale(i));
             end
            fprintf(funit,'\n');
            fprintf(1    ,'\n');
        end
        fclose(funit);    
    else
        error(sprintf('Could not open file named %s to record parameters.\n',fnameout));
    end
end
return
end

