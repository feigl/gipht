function [p1,cost1,trials,costs,kount,msig] = gridsearch1(objfun,bounds,OPTIONS,fitfun,DST,PST,TST)
%
%Version 1.1  Kurt Feigl 2012-NOV27

%pb=NaN;
cost1=NaN;
trials=NaN;
costs=NaN;
kount = 0;
msig = NaN;
talk = 1;

fprintf(1,'starting %s\n',mfilename);

%Check bounds to make sure they're ok
if max(bounds(:,1)>bounds(:,2))
    error('All the values in the first column of bounds must be less than those in the second.');
end


% initial estimate is center of bounds
%p0 = (bounds(:,1)+bounds(:,2))/2.0;
% initial estimate 
p0 = PST.p0;
% current best
pb = p0;
mparams = numel(p0);

% This is important
objfunhandle = str2func(char(objfun));
fitfunhandle = str2func(char(PST.fitfun));

% get initial cost
tstart=tic;
cost0=feval(objfunhandle,p0,fitfun,DST,PST,TST);
trun1 = toc(tstart);
fprintf(1,'Time for 1 evaluation is %.6f seconds\n',trun1);
cost1 = cost0;

% find indices of free parameters
ifree = find(bounds(:,2) > bounds(:,1));
fprintf(1,'Number of free parameters is %d\n',numel(ifree));
% number of slices in each parameter
nslice = 50;
%nslice = 10;

nruns = numel(ifree)*nslice;

fprintf(1,'Number of evaluations will be %d\n',nruns);
fprintf(1,'Time for %d evaluations will be %.1f seconds\n',nruns,nruns*trun1);

trials = zeros(nruns,mparams);
costs = nan(nruns,1);
fprintf(1,'kount  COST\n');
for ip=1:numel(ifree) 
    fprintf(1,'Varying parameter %d %s\n',ifree(ip),char(PST.names{ifree(ip)}));
    dp = (bounds(ifree(ip),2) - bounds(ifree(ip),1))/nslice;
    for is = 1:nslice
        kount = kount+1;
        
        % set up a trial set of parameters
        ptrial = pb;
        ptrial(ifree(ip)) = bounds(ifree(ip),1) + dp*(is-1);
        trials(kount,:) = ptrial;      
        
        % evaluate cost of trial
        ctrial=feval(objfunhandle,ptrial,fitfunhandle,DST,PST,TST);
        costs(kount) = ctrial;       
        
        if ctrial < cost1 
            fprintf(1,'%5d %12.4f\n',kount,ctrial);
            pb = ptrial;
            cost1 = ctrial;
        end
    end
end

if cost1 < cost0
    fprintf(1,'Found an estimate with cost %12.4e\n',cost1);
    p1 = pb;
else
    fprintf(1,'Did not find an estimate with lower cost.\n');
    fprintf(1,'Retaining initial estimate with cost %12.4e\n',cost0);
    p1 = p0;
end
fprintf(1,'leaving %s\n',mfilename);


return


