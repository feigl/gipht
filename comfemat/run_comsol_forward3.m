function modeledDisplacements = run_comsol_forward3(PST,DST,TST)
% PST parameter structure
% DST data structure
% fitting function
import com.comsol.model.util.*

narginchk(3,3);

verbose = 2;

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
        scale1 = PST.scale(ip);
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

%% run comsol solution and save it
fprintf(1,'%s: starting COMSOL solution with mph file named %s\n',mfilename,PST.fileNameMPH);
%fprintf(1,'%s: starting COMSOL solution...\n',mfilename);
trun0 = tic;
model.sol(PST.soltag).runAll;

mphsave(model);
fprintf(1,'%s: finished COMSOL solution in %10.1f seconds\n',mfilename,toc(trun0));


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

return
end

