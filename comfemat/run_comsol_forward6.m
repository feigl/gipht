function Lmod = run_comsol_forward6(PST,DST,TST)
% PST parameter structure
% DST data structure
% TST temporary structure
% fitting function
import com.comsol.model.util.*

narginchk(3,3);

verbose = 2;
Lmod = nan(size(DST.obs));

%% load the mph file
model = mphload(PST.fileNameMPH);
%https://www.comsol.com/forum/thread/128742/avoid-mphload-every-iteration-in-matlab-livelink
% Inside you objective function you can now do of one these:
% 1)
% import com.comsol.model.util.*
% model = ModelUtil.model('MyModel');
% gives you access to the MyModel model on the server using the variable name model in Matlab

% count parameters
if numel(PST.p1) ~= numel(PST.name)
    error('miscount of parameters');
else
    mparams = numel(PST.p1);   
end

%% set parameters
nreset = 0;
for ip=1:mparams  % loop over indices to comsol parameters
    current1 = PST.p1(ip);  % current value 
    if isfinite(current1) == 1
        db = PST.ub(ip)-PST.lb(ip); % range of bounds in gipht parameter
        %scale1 = PST.scale(ip);
        if  db > 0
            %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
            pname1 = sprintf('%s',PST.name{ip});
            unit1 = char(model.param.evaluateUnit(pname1));
            
            if strcmp(unit1,PST.unit(ip)) == false
                fprintf(1,'Unit in MPH file (%s) does not match unit in PST structure (%s) for parameter named %s\n'...
                    ,unit1,PST.unit(ip),pname1);
            end
            desc1 = char(model.param.descr(pname1));
            if verbose >= 2
                %fprintf(1,'    Setting %d %-8s to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),PST.description{ip});
                fprintf(1,'    Setting %d %-8s to %12.4e [%s] %s\n',ip,pname1,current1,unit1,desc1);
            end
            %model.param.set(pname1,sprintf('%e [%s]',param1,unit1),desc1);
            % rescale from dimensionless back to units in comsol
            model.param.set(pname1,sprintf('%e [%s]',current1,unit1),desc1);
            nreset = nreset + 1;
        end        
    end
end

% get some information
model_struct = mphmodel(model);
studys = model_struct.study;
nstudy=numel(studys);

%model=mphload(fileNameMPH)


%% run the solution
trun0 = tic;
fprintf(1,'%s: starting COMSOL solution with mph file named %s\n',mfilename,PST.fileNameMPH);
for i=1:nstudy
    study1 = studys{i};
    model.study(study1).run;
end
mphsave(model);
fprintf(1,'%s: finished COMSOL solution in %10.1f seconds\n',mfilename,toc(trun0));


%% calculate modeled values at locations of data points
ndata = numel(DST.obs);

% find unique locations by distance from orgin
rdist=sqrt(DST.x.^2 + DST.y.^2 + DST.z.^2);
[runique,iunique] = unique(rdist);
nObsPoints = numel(runique);
XYZobspts=zeros(3,nObsPoints);

XYZobspts(1,:) = reshape(DST.x(iunique),1,nObsPoints);
XYZobspts(2,:) = reshape(DST.y(iunique),1,nObsPoints);
XYZobspts(3,:) = reshape(DST.z(iunique),1,nObsPoints);

% times when data occur in seconds elapsed from start of model
nEpochs = numel(DST.t)
%tEpochs = reshape(DST.t - PST.t0,1,nEpochs);
tEpochs = reshape(DST.t,1,nEpochs);

%% extract modeled displacements by interpolating coordinates and time
Umod = mphinterp(model,{'u2'},'coord',XYZobspts,'t',tEpochs,'dataset',PST.dset); % easting component
Vmod = mphinterp(model,{'v2'},'coord',XYZobspts,'t',tEpochs,'dataset',PST.dset); % northing component
Wmod = mphinterp(model,{'w2'},'coord',XYZobspts,'t',tEpochs,'dataset',PST.dset); % vertical component


[nur,nuc] = size(Umod);
if nur ~= nEpochs || nuc ~= nObsPoints 
    nur
    nuc
    nEpochs
    nObsPoints
    warning('arrays are the wrong size');
end

% make displacement relative to first observed epoch 
for i=1:nObsPoints
    Umod(:,i) = Umod(:,i) - Umod(1,i);
    Vmod(:,i) = Vmod(:,i) - Vmod(1,i);
    Wmod(:,i) = Wmod(:,i) - Wmod(1,i);
end


% compute displacement along line of sight
Lmod = Umod.*DST.pointingX + Vmod.*DST.pointingY + Wmod.*DST.pointingZ;





return
end

