%function modeledDisplacements = run_comsol_forward5(PST,DST,TST)
% PST parameter structure
% DST data structure
% TST temporary structure
% 20200627 Kurt Feigl - handle-time dependent solution
import com.comsol.model.util.*

%clear all
% narginchk(3,3);

%% Read GPS velocity data
%file='/Users/hlemevel/Documents/MATLAB/3 dMODELS (2019)/1 GPS/MATLAB/sill/LDM_2019_2020_GPS/LdM_Dmodel_GPSvel20192020.txt';
file='LdM_Dmodel_GPSvel20192020.txt';
format='%s %f %f %f %f %f %f %f %f %f';
[sitesObs,xutm,yutm,alt,E,dE,N,dN,U,dU] = textread(file,format,'delimiter',' ','headerlines',3);

%% This list MUST match those in the model
DST.sites = {'MAU2','PUEL','COLO','LDMP','NIE2'}
nsites = numel(DST.sites)
ndata = 3 * nsites;
DST.obs = nan(ndata,1);
DST.sig = nan(ndata,1);
DST.x = nan(ndata,1);
DST.y = nan(ndata,1);

%% get observed values and measurement uncertainties
idata=0;
for i=1:nsites
    isite=find(contains(sitesObs,upper(DST.sites{i})));
    if numel(isite) == 1
        idata=idata+1; DST.x(idata) = xutm(isite); DST.y(idata) = yutm(isite); DST.obs(idata) = E(isite); DST.sig(idata) = dE(isite);
        idata=idata+1; DST.x(idata) = xutm(isite); DST.y(idata) = yutm(isite); DST.obs(idata) = N(isite); DST.sig(idata) = dN(isite);
        idata=idata+1; DST.x(idata) = xutm(isite); DST.y(idata) = yutm(isite); DST.obs(idata) = U(isite); DST.sig(idata) = dU(isite);
    else
        error('sites not in correct order');
    end
end
if idata ~= ndata
    error('miscount reading data file');
end


verbose = 2;

%% load the mph file
%PST.fileNameMPH = 'LdM_3DFSI_spheroid_Tvisco_P_forward2020JUN25save3.mph'
%PST.fileNameMPH = '/Users/feigl/BoxSync/LDMcomsol/LdM_3DFSI_spheroid_Tvisco_FORWARD2020MAY06nosweep.mph'
%PST.fileNameMPH = 'LdM_3D_ellipsoid_FORWARD_shape_2020JUNE5.mph'
% PST.soltag='sol4'
% PST.dset = 'dset9'
% PST.solnum = 2

%model = mphload(PST.fileNameMPH);
%https://www.comsol.com/forum/thread/128742/avoid-mphload-every-iteration-in-matlab-livelink
% Inside you objective function you can now do of one these:
% 1)
% import com.comsol.model.util.*
% model = ModelUtil.model('MyModel');
% gives you access to the MyModel model on the server using the variable name model in Matlab

% %% get information
% %info1 = mphsolinfo(model,'soltag','sol1')
% info0 = mphsolutioninfo(model)
% for i=1:numel(info0.solutions)
%     info1 = mphsolinfo(model,'soltag',info0.solutions{i})
% end
% 
% %%
% if verbose == 1
%     info0
%     info1
% end

%field1 = mphgetfield(model.field('field3'))

% field1 = 
% 
%   struct with fields:
% 
%     field: 'u2'
%     shape: 'shape5, shape6, shape7'
%      geom: 'geom1'
 

% Comsol evaluation solution_epochs in seconds
%solution_epochs = info1.solvals;

% verbose = 0;
%Tparams = get_comsol_parameters2(PST.fileNameMPH,verbose);
Tparams = get_comsol_parameters2(model,verbose);
%% read file of displacements in UNR 'tenv3' format
ip=find(strcmp(Tparams.name,'t0'));DST.t0 = Tparams.value(ip) % initial time in seconds


% % count parameters
% if numel(PST.p1) ~= numel(PST.name)
%     error('miscount of parameters');
% else
%     mparams = numel(PST.p1);   
% end

% %% set parameters
% nreset = 0;
% for ip=1:mparams  % loop over indices to comsol parameters
%     current1 = PST.p1(ip);  % current value 
%     if isfinite(current1) == 1
%         db = PST.ub(ip)-PST.lb(ip); % range of bounds in gipht parameter
%         scale1 = PST.scale(ip);
%         if  db > 0
%             %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
%             pname1 = sprintf('%s',PST.name{ip});
%             unit1 = char(model.param.evaluateUnit(pname1));
%             
%             if strcmp(unit1,PST.unit(ip)) == false
%                 fprintf(1,'Unit in MPH file (%s) does not match unit in PST structure (%s) for parameter named %s\n'...
%                     ,unit1,PST.unit(ip),pname1);
%             end
%             desc1 = char(model.param.descr(pname1));
%             if verbose >= 2
%                 %fprintf(1,'    Setting %d %-8s to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),PST.description{ip});
%                 fprintf(1,'    Setting %d %-8s to %12.4e [%s] %s\n',ip,pname1,current1,unit1,desc1);
%             end
%             %model.param.set(pname1,sprintf('%e [%s]',param1,unit1),desc1);
%             % rescale from dimensionless back to units in comsol
%             model.param.set(pname1,sprintf('%e [%s]',current1,unit1),desc1);
%             nreset = nreset + 1;
%         end        
%     end
% end
% 
% %% run comsol solution and save it
% fprintf(1,'%s: starting COMSOL solution with mph file named %s\n',mfilename,PST.fileNameMPH);
% %fprintf(1,'%s: starting COMSOL solution...\n',mfilename);
% trun0 = tic;
% model.sol(PST.soltag).runAll;
% 
% mphsave(model);
% fprintf(1,'%s: finished COMSOL solution in %10.1f seconds\n',mfilename,toc(trun0));



% Tgps = readtable('GARL.tenv3','filetype','text');
% Tgps(1:10,:)
% DST.sites = {'GARL'};
% DST.x = Tgps.x_e0_m_    + Tgps.x__east_m_;  % UTM easting in meters
% DST.y = Tgps.x____n0_m_ + Tgps.x_north_m_;  % UTM northing in meters
% DST.z = Tgps.u0_m_      + Tgps.x____up_m_;  % UTM elevation in meters (assumed above WGS84 ellipsoid)
% DST.t = (Tgps.yyyy_yyyy*3600.*24.*365.25) - DST.t0;   % elapsed time in seconds

% mphinterp(model,{'u2'},'coord',[DST.x';DST.y';DST.z'],'t',DST.t,'solnum','all')
% Error using mphinterp
% Unable to perform assignment because the size of the left side is 0-by-8499 and the size of the right side is 1-by-8499.
%  
%mphinterp(model,{'u2'},'coord',[DST.x';DST.y';DST.z'],'t',DST.t,'outersolnum','all')
% Error using mphinterp
% Unable to perform assignment because the size of the left side is 0-by-8499 and the size of the right side is 1-by-8499.
% Error using mphinterp
% Java exception occurred:
% Exception:
% 	com.comsol.util.exceptions.FlException: Variable refers to an object which is no longer in the model
% Messages:
% 	Variable refers to an object which is no longer in the model.


info0 = mphsolutioninfo(model)

%% calculate modeled values at locations of data points
%expr = mphgetexpressions(model.param);
sites = DST.sites;
nsites = numel(sites);
ndata = 3 * nsites;
modeledDisplacements = nan(ndata,1);
idata=0;

%% calculate modeled values at locations of data points
expr = mphgetexpressions(model.param);
for isite = 1:nsites
    ix=find(contains(expr(:,1),lower(sprintf('x%4s',sites{isite}))));
    iy=find(contains(expr(:,1),lower(sprintf('y%4s',sites{isite}))));
    epochs = DST.t0 + (0:20)*24.*3600.*365.25; % in absolute seconds after time 0
    nepochs = numel(epochs)
    coords = [expr{ix,4};expr{iy,4};0];
    [ncoords,ncols] = size(coords)
    %coords=repmat([expr{ix,4};expr{iy,4};0],1,nepochs);
    
    % %% get information
    % %info1 = mphsolinfo(model,'soltag','sol1')
    %
    % %%
    % if verbose == 1
    %     info0
    %     info1
    % end
    
    for isol=1:numel(info0.solutions)
    %for isol=5
        info1 = mphsolinfo(model,'soltag',info0.solutions{isol});
        dsets = info1.dataset;
        if isempty(dsets) == true
            ndsets = 0;
        else
            if iscell(dsets) == true
                ndsets = numel(dsets);
            else
                ndsets = 1;
            end
        end
        for idset = 1:ndsets
            if ndsets > 1
                dset1 = dsets{idset};
            else
                dset1 = dsets;
            end
            
            clear U;
            try
                fprintf(1,'Trying solnum %d data set %s\n',isol,dset1);
                %U = mphinterp(model,{'u2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1,'t',DST.t0 + ([0:20]*24.*3600.*365.25),'solnum',isol); %east
                %U = mphinterp(model,{'u2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1,'t',DST.t0 + ([0:20]*24.*3600.*365.25),'outersolnum',isol); %east
                %U = mphinterp(model,{'u2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1,'t',DST.t0 + ([0:20]*24.*3600.*365.25)); %east
                %size(mphinterp(model,{'u2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1,'t',DST.t0 + ([0:20]*24.*3600.*365.25))) %east
                %U = mphinterp(model,'u2','coord',coords,'dataset',dset1,'t',epochs); %east displacement
                X = mphinterp(model,'x','coord',coords,'dataset',dset1,'t',epochs); %easting coordinate
%                 if numel(U) > 0
%                     whos U
%                     fprintf(1,'Success with solnum %d data set %s\n',isol,dset1);
%                 end
                
            catch ME
                ME.message
                warning('mphinterp failed');
            end
            
            %         V = mphinterp(model,{'v2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1); % north component
            %         W = mphinterp(model,{'w2'},'coord',[expr{ix,4};expr{iy,4};0],'dataset',dset1); % vertical component
            %
            
            %     U = mphinterp(model,{'u2'},'coord',[DST.x';DST.y';DST.z'],'t',DST.t,'dataset',PST.dset,'outersolnum',PST.solnum); % eastward component of displacement in meters
            %     V = mphinterp(model,{'v2'},'coord',[DST.x';DST.y';DST.z'],'t',DST.t,'dataset',PST.dset,'solnum',PST.solnum); % north component of displacement in meters
            %     W = mphinterp(model,{'w2'},'coord',[DST.x';DST.y';DST.z'],'t',DST.t,'dataset',PST.dset,'solnum',PST.solnum); % vertical component of displacement in meters
            %     idata=idata+1;modeledDisplacements(idata) = U;
            %     idata=idata+1;modeledDisplacements(idata) = V;
            %     idata=idata+1;modeledDisplacements(idata) = W;
            
        end
        
    end
    
    
%     if idata ~= ndata
%         error('miscount')
%     end
    
    return
end

