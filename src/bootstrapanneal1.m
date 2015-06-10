function [p1,psigma,cost1,pboot]=bootstrapanneal1(objfun,bounds,p0,options,varargin)
%function [p1,psigma,cost1]=bootstrapanneal1(objfun,bounds,p0,options,varargin)
%
% perform bootstrap analysis to find uncertainties of parameter estimates
% This version uses simulated annealing.
%
% Efron, B., and R. Tibshirani (1986), 
% Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy, 
% Statistical
% Science, 1, 54-77.  
% inputs:
%        objfun      == name of objective function
%        bounds   == matrix with 2 rows, giving lower and upper bounds of
%                    each parameter
%        options  == vector of options
%
% Version 1.1  Kurt Feigl 2012-JAN-06


%Check argument syntax

if nargin<2
    help(mfilename);
    error('Incorrect arguments');
end


if size(bounds,2)~=2
    error('Second argument must be an nx2 matrix of parameter bounds, where n is the number of parameters.');
end

p1 = nan(size(p0));
psigma = nan(size(p0));
pmean = nan(size(p0));
cost1 = NaN;


%Check options
if nargin<4
    options=[];
end
if isempty(options)
    matrix=0;
    newton=0;
    talk=1;
    parallel=0;
else
    matrix=options(5);
    newton=options(6);
    if numel(options) >= 7
        talk=options(7);
    end
    if numel(options) >= 8
        if options(8) > 1
            fprintf(1,'Checking that matlabpool is open with %d processors....\n',options(8));
            %matlabpool('open',options(8))
            if matlabpool('size') == options(8)
                fprintf(1,'Success.\n');
            else
                error(sprintf('Distproc failure.\n'));
            end
        end
    end
end



DST = varargin{2};
PST = varargin{3};
TST = varargin{4};

objfunhandle = str2func(char(objfun));
fitfunhandle = str2func(char(PST.fitfun));

% evaluate cost of initial estimate
%cost0=feval(objfun,p0,varargin{:});
cost0 = feval(objfunhandle,p0,fitfunhandle,DST,PST,TST);

if talk
    fprintf('\nStarting Jacknife using Simulated annealing.\n\n')
end

DST0 = DST;

% number of parameters
mparam = numel(p0);

% Number of bootstrap resamples
nboot = 100; % "There is little improvement past B = 100"
% nboot = 3; % test quickly
pboot = nan(nboot,mparam);
costsj = nan(nboot,1);

% perturb initial estimate slightly
%p0j = p0 + 0.1 * abs(bounds(:,2)-bounds(:,1));

fprintf(1,'i fval nanmean(DST.phaobs-DST.phamod) nanstd(DST.phaobs-DST.phamod)\n');
for i = 1:nboot
    DST = bootstrap_dst(DST0);
    
    % re-make TST structure
    [dummy,TST] = feval(fitfunhandle,DST,PST);
    
    [pj,costj,trials,energy,count] = anneal5(objfunhandle,bounds,options,fitfunhandle,DST,PST,TST);

    
    if sum(isfinite(pj)) == mparam
        pboot(i,:)=pj;
        costsj(i)= costj;
        fprintf(1,'%5d %20.10E %20.10E %20.10E\n',i,costj,nanmean(DST.phaobs-DST.phamod),nanstd(DST.phaobs-DST.phamod));
    end
end

% Jackknife standard deviation from equation (4)
fprintf(1,'j,p1(j),psigma(j)\n');
mfree = 0;
for j=1:mparam
    if abs(bounds(j,2)-bounds(j,1)) > 0
        mfree = mfree+1;
        psigma(j) = nanstd(pboot(:,j));
        pmean(j) = nanmean(pboot(:,j));
        fprintf(1,'%5d %12.4E +/- %12.4E\n',j,pmean(j),psigma(j));
    end
end


% retrieve best values
[cost1,ibest] = nanmin(costsj);
p1 = pboot(:,ibest);
return


