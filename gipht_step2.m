% GIPHT_step2
% Set bounds on parameters
% Run simulated annealing
%
% 2009-APR-13 slight modifications to case with anneal == 2
% 2009-JUN-18 integer phase
% 2009-OCT-21 handle individual epochs or wild card (epoch_00)

% For stand alone executable version, replace feval with feval2?
fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));


%clear all
clearvars

if fexist('gipht.mat') == 1
    load('gipht.mat');
end

% Real and imaginary parts of phase
%xd = cos(phao*2*pi); % argument in radians
%yd = sin(phao*2*pi); % argument in radians
%xd = phao;  % phase in DN [-128, 127] 256 DN = 1 cycle
%xd = phao; % phase in radians [-pi, +pi] 2010 JAN 11
%xd = DST.phaobs;
%whos xd
%yd = NaN * ones(size(xd)); % for historical reasons
%yd = zeros(size(xd),'int8'); % for historical reasons
%yd = zeros(size(xd)); % double 2010 JAN 11

ndata = numel(DST.phaobs)
ndata1 = DST.i(ndata)
if ndata1 ~= ndata 
    error(sprintf('ERROR: data miscount ndata1 (%d) does not equal ndata (%d)!\n',ndata1,ndata));
end


% get names of parameters and scaling factors
[pnames,pscl] = feval(fitfun,me,datafilename);
params_per_epoch = 11;
disp('Number of parameters')
mparam = numel(pnames)
p0=zeros(mparam,1);

% bounds: fix all parameters
lb=zeros(mparam,1);
ub=zeros(mparam,1);

% read initial and bounding values of parameters from file if extant
fidparin = fopen(fnparin,'r');
if fidparin > 0
   fprintf(1,'Opened file named %s to read initial estimates of parameters.\n',fnparin);
   for i = 1:mparam
       parse_err=0;
       k=0;
       frewind(fidparin);
       aline = fgetl(fidparin);
       while 1
           aline = fgetl(fidparin);
           % 2011-NOV-30 if aline < 0
           % 2011-NOV-30 allow blank lines
           if numel(aline) < 40
               parse_err = -1;
               break;
           end
           
           scanned=textscan(aline,'%s%s%s','Delimiter',' ','MultipleDelimsAsOne',1);
           pnamein = scanned{1};
           pnamein = char(scanned{1}); % 2009-OCT-27 convert cell to string
           if i <= params_per_epoch * me
               pnamein2 = strrep(pnamein,   '_epoch_000',sprintf('_epoch_%03d',loopfun(i,me))); pnamein = pnamein2;
               pnamein2 = strrep(pnamein,   '_epoch_00_',sprintf('_epoch_%03d',loopfun(i,me))); pnamein = pnamein2;
           end
           %fprintf(1,'Parsing parameter indexed %4d %s %s\n',i,char(pnames{i}),pnamein);
          
           if regexpindex(pnames,pnamein) == i   % 2009-OCT-21
               %fprintf(1,'Found parameter indexed %4d %s %s\n',i,char(pnames{i}),pnamein);
               valinit  = str2double(char(scanned{2}));
               if isfinite(valinit) == 1
                   p0(i) = valinit;
               else
                   p0(i) = NaN;
                   parse_err = parse_err + 1;
               end
               plusminus=str2double(char(scanned{3}));
               if isfinite(plusminus) == 1
                   if isfinite(p0(i)) == 1
                       lb(i) = p0(i) - plusminus;
                       ub(i) = p0(i) + plusminus;
                   else
                       parse_err = parse_err + 1;
                   end
               else
                   parse_err = parse_err + 1;
               end
               if parse_err ~= 0
                   fprintf(1,'ERROR: garbled parameter. Expected:\n%s p0_value plus_minus\nFound:\n%s %s %s\n'...
                       ,char(pnames{i}),char(pnamein),char(scanned{2}),char(scanned{3}));
                   fprintf(1,'1LB P0 UB for %3d %s %12.6g %12.6g %12.6g\n',i,pnames{i},lb(i),p0(i),ub(i));
                   error(sprintf('parse_err = %d',parse_err));
               end
           end
       end % while loop over lines of file
   end % loop over parameters
   fclose(fidparin);
   fprintf(1,'\n Done reading parameters\n');

%% TODO read table and make into a structure
%    GINT = readtable('demoF14.gin','filetype','text','ReadVariableNames',true,'HeaderLines',0)
% GINT = 
%     model_name    ID             parameter_name               initial     plusminus
%     __________    __    _________________________________    _________    _________
%     'EPOCH'       0     'Bperp___@_epoch_000_in_m_______'            0       0     
%     'EPOCH'       0     'E_grad__@_epoch_000_dimless____'            0       0     
%     'EPOCH'       0     'N_grad__@_epoch_000_dimless____'            0       0     
%     'EPOCH'       0     'U_grad__@_epoch_000_dimless____'            0       0     
%     'OKADA'       1     'Centroid_Easting_in_m__________'     5.09e+05     300     
%     'OKADA'       1     'Centroid_Northing_in_m_________'    3.802e+06     300     
%     'OKADA'       1     'Centroid_Depth_in_m____________'         2000     500     
%     'OKADA'       1     'Width_in_m_____________________'         2340     300     
%     'OKADA'       1     'Length_in_m____________________'         2780     300     
%     'OKADA'       1     'Dip_in_deg____________________'            50      10     
%     'OKADA'       1     'Strike_CW_from_N_in_deg________'          283       5     
%     'OKADA'       1     'Coplanar_slip_in_m_____________'         0.52    0.05     
%     'OKADA'       1     'Rake_in_deg_CCW________________'           70      10     
%     'OKADA'       1     'Tensile_Opening_in_m___________'            0       0     
%     'OKADA'       1     'Poisson_Ratio_dimless__________'         0.25       0     
%     'OKADA'       1     'Shear_Modulus_in_Pa____________'        3e+10       0     
%     'TFUNC'       1     'Reference_Epoch_in_years_______'       1992.9       0     
% PST = table2struct(GINT,'ToScalar',true)
% PST = 
%         model_name: {17x1 cell}
%                 ID: [17x1 double]
%     parameter_name: {17x1 cell}
%            initial: [17x1 double]
%          plusminus: [17x1 double]
% PST.parameter_name(1)
else
    warning(sprintf('Could not open input file: %s. Using zeros\n',txtinname));
end % reading parameters from file

% % Additive constant depends on observable
%20170802 if ismember(pselect,[7,9]) == 1  % observable is gradient
switch idatatype1
    case 0  % observable is phase
        p0(get_parameter_index('Offset',pnames)) =  0.0; % cycles
        lb(get_parameter_index('Offset',pnames)) = -0.5; % cycles
        ub(get_parameter_index('Offset',pnames)) = +0.5; % cycles
        ifast = 0;
    case 2 % observable is range - bounds to be defined later
%         p0(get_parameter_index('Offset',pnames)) =  nan;
%         lb(get_parameter_index('Offset',pnames)) =  nan;
%         ub(get_parameter_index('Offset',pnames)) =  nan;
%       allow for constant value 
        p0(get_parameter_index('Offset',pnames)) =  0.0; % meters?
        lb(get_parameter_index('Offset',pnames)) = +nanmean([nanmin(DST.phaobs),nanmax(DST.phaobs)]);
        ub(get_parameter_index('Offset',pnames)) = -nanmean([nanmin(DST.phaobs),nanmax(DST.phaobs)]);
        ifast = 0;
    case -1 % observable is gradient
        p0(get_parameter_index('Offset',pnames)) =  0.0; % cycles
        lb(get_parameter_index('Offset',pnames)) =  0.0; % cycles
        ub(get_parameter_index('Offset',pnames)) =  0.0; % cycles
        ifast = 100;
    otherwise
        error(sprintf('Unknown idatatype1 %d\n',idatatype1));
end
%

% set initial values of epoch-wise parameters
for k=1:me  %    % loop over epochs
    % index of master
    %kmast = find(DD(k,:) == -1);
    % index of pair
    kp = find(DD(:,k) == -1); % epoch is master
    if numel(kp) == 0
        kp = find(DD(:,k) == +1); % epoch is slave
    end
    if numel(kp) == 0
        error('pair index not found!');
    end
    % index to pixels
    i1 = ippix1(kp);
    if kp < np
        i2 = ippix1(kp+1);
    else
        i2 = ndata;
    end
    
    [is,jm] = find(trees == k);
    % is is index to trees
    % jm is index of member within trees
    
    %pindex=             get_parameter_index(sprintf('time_fn_@_epoch_%03d_____________',k),pnames);
    pindex = [];%                                     E_grad__@_epoch_001_dimless_____
    % 2012-JAN-10 next line is a bug!!
    %pindex=union(pindex,get_parameter_index(sprintf('Offset__@_epoch_%03d_____________',k),pnames));
    %pindex=union(pindex,get_parameter_index(sprintf('Offset__@_epoch_%03d_in_cycles___',k),pnames));
    % 20171213 allow offset in other units
    pindex=union(pindex,get_parameter_index(sprintf('Offset__@_epoch_%03d______________',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('E_grad__@_epoch_%03d_dimless_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('N_grad__@_epoch_%03d_dimless_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('U_grad__@_epoch_%03d_dimless_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitHoriz_@_epoch_%03d_in_m_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitAlong_@_epoch_%03d_in_m_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitVerti_@_epoch_%03d_in_m_____',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitVelH__@_epoch_%03d_m_per_s__',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitVelA__@_epoch_%03d_m_per_s__',k),pnames));
    pindex=union(pindex,get_parameter_index(sprintf('OrbitVelV__@_epoch_%03d_m_per_s__',k),pnames));
    
    for ii=1:numel(pindex)
        if jm == 1   % epoch is first in its tree
            %if abs(p0(pindex(ii))) >= 0.0 && isfinite(p0(pindex(ii))) == 1 % initial value is already set
            if isfinite(p0(pindex(ii))) == 1 % initial value is already set
                lb(pindex(ii)) = p0(pindex(ii));                % update lower bound
                ub(pindex(ii)) = p0(pindex(ii));                % update upper bound
                fprintf(1,'Fixing bounds to initial estimate on %s %16.7g %16.7g %16.7g\n'...
                    ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
            else % initial value is not already set
                p0(pindex(ii))=  0.0;
                lb(pindex(ii)) = 0.0;                % update lower bound
                ub(pindex(ii)) = 0.0;                % update upper bound
                fprintf(1,'Zeroing  initial value and bounds on %s %16.7g %16.7g %16.7g\n'...
                    ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
            end
        else % epoch is NOT first in its tree
% %             % 20170802 add this case
%             if isfinite(p0(pindex(ii))) == 1 % initial value is already set
%                 lb(pindex(ii)) = p0(pindex(ii));                % update lower bound
%                 ub(pindex(ii)) = p0(pindex(ii));                % update upper bound
%                 fprintf(1,'Fixing bounds to initial estimate on %s %16.7g %16.7g %16.7g\n'...
%                     ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
%             else               
                lb(pindex(ii)) = p0(pindex(ii)) + lb(pindex(ii));                % update lower bound
                ub(pindex(ii)) = p0(pindex(ii)) + ub(pindex(ii));                % update upper bound
                fprintf(1,'Updating                       bounds on %s %16.7g %16.7g %16.7g\n'...
                    ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
 %           end
        end
    end
    
    % should not estimate gradient parameters
    if ismember(pselect,[7,9]) == 1  % observable is gradient
        pindex=get_parameter_index(sprintf('E_grad__@_epoch_%03d_dimless_____',k),pnames);
        for ii=1:numel(pindex)
            if isfinite(p0(pindex(ii))) == 1 % initial value is already set
                lb(pindex(ii)) = p0(pindex(ii));                % update lower bound
                ub(pindex(ii)) = p0(pindex(ii));                % update upper bound
                fprintf(1,'Fixing bounds to initial estimate on %s %16.7g %16.7g %16.7g\n'...
                    ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
            else % initial value is not already set
                p0(pindex(ii))=  0.0;
                lb(pindex(ii)) = 0.0;                % update lower bound
                ub(pindex(ii)) = 0.0;                % update upper bound
                fprintf(1,'Zeroing  initial value and bounds on %s %16.7g %16.7g %16.7g\n'...
                    ,pnames{pindex(ii)},p0(pindex(ii)),lb(pindex(ii)),ub(pindex(ii)));
            end
        end
    end
end

%% handle time function
pindex=get_parameter_index('time_fn',pnames);
if numel(pindex) == numel(tepochs)
%     decyears = dyear(year(tepochs),month(tepochs),day(tepochs));
%     p0(pindex) = decyears;
%     lb(pindex) = decyears;
%     ub(pindex) = decyears;
    decyears = dyear(year(tepochs),month(tepochs),day(tepochs));
    p0(pindex) = decyears;
    lb(pindex) = decyears;
    ub(pindex) = decyears;
else
    error('numel(pindex) NE numel(tepochs)');
end



% count number of degrees of freedom
if (numel(ub) == numel(lb)) && (numel(p0) == numel(ub))
    mparam = numel(p0);
else
    error 'ERROR: problem with number of parameters'
end

% display parameter values
fprintf(1         ,'I          Name     P0(pre-fit)              P1(post-fit)     Adjust.   Sigma[%1d] Significance PlusMinusBound\n',istatcode);

for i = 1:mparam
%     if numel(findstr(pnames{i},'time_fn')) > 0 || ( abs(p0(i)) > 0.1 && abs(p0(i)) < 10000)
%         fmtstr = '%3d %#28s %10.3f %10.3f %10.3f %8s\n';
%     else
%         fmtstr = '%3d %#28s %10.2e %10.2e %10.2e %8s\n';
%     end
    
    if abs(ub(i)-lb(i)) > 1e-6 * p0(i)
        %bstat = 'free ';
        pflag{i} = 'E#';
    else
        %bstat = 'fixed';
        pflag{i} = 'F#';
    end
    psig(i) = NaN;
    adj = 0;
    sadj    = NaN;
    p1(i)   = NaN;
    
    %fprintf(1,fmtstr,i,pnames{i},p0(i),lb(i),ub(i),bstat);
    
    outfmt = getfmt(p0(i),pnames{i});
    fprintf(1        ,outfmt,pflag{i},i,pnames{i} ,p0(i),p1(i),adj,psig(i),sadj,(ub(i)-lb(i))/2.0);
    
    % check values here
    if findstr(lower(pnames{i}),'easting') > 0
        if (p0(i)/1.0e0 < xsubmin || p0(i)/1.0e0 > xsubmax) && abs(p0(i)) > 1e-6
            fprintf('WARNING: initial estimate of easting coordinate outside subregion.\n')
        end
        if (lb(i)/1.0e0 < xsubmin || lb(i)/1.0e0 > xsubmax) && abs(ub(i)-lb(i)) > 1e-6
            fprintf('WARNING: lower bound on easting coordinate outside subregion.\n')
        end
        if (ub(i)/1.0e0 < xsubmin || ub(i)/1.0e0 > xsubmax) && abs(ub(i)-lb(i)) > 1e-6
            fprintf('WARNING: upper bound on easting coordinate outside subregion.\n')
        end
    end
    if findstr(lower(pnames{i}),'northing') > 0
        if (p0(i)/1.0e0 < ysubmin || p0(i)/1.0e0 > ysubmax) && abs(p0(i)) > 1e-6
            fprintf('WARNING: initial estimate of northing coordinate outside subregion.\n')
        end
        if (lb(i)/1.0e0 < ysubmin || lb(i)/1.0e0 > ysubmax) && abs(ub(i)-lb(i)) > 1e-6
            fprintf('WARNING: lower bound on northing coordinate outside subregion.\n')
        end
        if (ub(i)/1.0e0 < ysubmin || ub(i)/1.0e0 > ysubmax) && abs(ub(i)-lb(i)) > 1e-6
            fprintf('WARNING: upper bound on northing coordinate outside subregion.\n')
        end
    end
    if findstr(lower(pnames{i}),'U_grad') > 0
        if pselect == 3 && abs(ub(i)-lb(i)) > 1e-6
            fprintf('WARNING: parameter %s cannot be estimated because DEM not read for PSELECT=3\n', pnames{i});
            fprintf('WARNING: Setting bounds on above parameter to be %.g4 +/- 0.0\n', pnames{i},p0(i));
            lb(i)=p0(i);
            ub(i)=p0(i);
        end
    end
end


mfree = 0;
for i=1:mparam
    if ub(i)-lb(i) > 1.0e-5 * abs(p0(i))
        mfree = mfree+1;
    elseif ub(i)<lb(i)
        fprintf(1,'WARNING: Bad bounds with parameter %d %s\n%e %e %e\n'...
            ,i,pnames{i},p0(i),lb(i),ub(i));
        fprintf('WARNING: parameter %s cannot be estimated because upper bound is less than lower bound.\n', pnames{i});
        fprintf('WARNING: Setting bounds on above parameter to be %.g4 +/- 0.0\n', pnames{i},p0(i));
        lb(i)=p0(i);
        ub(i)=p0(i);
    end
end
fprintf(1,'Number of free parameters MFREE = %d\n',mfree);
if mfree == 0 && ianneal ~= 0
    warning('No parameters are adjustable. Skipping optimization. Setting anneal to zero.\n');
    ianneal = 0;
end
    
%% bounds matrix
% must be an nx2 matrix of bounds on parameter
% where n is the number of parameters.
bounds(:,1) = lb;
bounds(:,2) = ub;
%size(bounds)

%% build parameter structure PST00 for null model
psig = NaN*ones(size(p0));
for i=1:numel(p0)
    pflags{i} = 'N#';
end
PST00 = build_pst(fitfun,mparam,zeros(size(p0)),zeros(size(p0)),zeros(size(p0))...
    ,pnames,bounds,datafilename,pscl,pflags,timefun);


%% build parameter structure PST0 for initial estimate
p1 = p0;
PST0 = build_pst(fitfun,mparam,p0,p1,psig,pnames,bounds,datafilename,pscl,pflags,timefun);
ierr = check_struct(PST0);
%%ierr2 = write_pst(PST0,fnamepstin);


% call fitting function first time to initialize
% temporary storage structure TST
ierr = check_struct(DST);
clear rng
clear TST
objfun
fprintf(1,'Calling fitting function for first time to initialize.\n');
[rng,TST] = feval(fitfun,DST,PST0);
% 20150519 Try to avoid using feval
%fhandle = @fitfun;
% fhandle = @funfit28
% [rng,TST] = fhandle(DST,PST);
DST.phamod = rng;

% evaluate function at initial estimate
clear mdl0;
fprintf(1,'Calling fitting function for second time\n');
mdl0 = feval(fitfun,DST,PST0,TST);
% 20150519 Try to avoid using feval
% fhandle = @fitfun;
% mdl0 = fhandle(DST,PST0,TST);
DST.phamod = colvec(mdl0);
% 
% % make a 3-D plot using initial estimate
% nf=nf+1;
% h(nf) = plot_obsmod3d(DST);
% feval(printfun,sprintf('%s_OBSMOD3D',runname));
% 
% % make plot of observed versus initial model
% nf=nf+1;h(nf)=figure;subplot(2,1,1);
% plot(rwrapm(DST.phaobs)/2.0/pi,'ro');hold on;
% plot(rwrapm(DST.phamod)/2.0/pi,'k-');
% title('initial estimate');
% ylabel('phase (cycles)');
% legend('observed','modeled');
% subplot(2,1,2);
% plot(rwrapm(DST.phaobs-DST.phamod)/2.0/pi,'bo-');
% xlabel('pixel index');ylabel('phase (cycle)');legend('residual');
% feval(printfun,sprintf('%s_Taylor1Check',runname));

% % Approximate Fitting function by a surrogate 
% % generate partials for first order Taylor series approximation
% if surrogate == 1
%     imode = 1;
%     maxiter = 1;
%     [TSTP, PST_dummy, trials] = anneal_with_surrogate(objfun, DST, PST, TST, maxiter ...
%         , verbose, imode, nf, runname, printfun);
%     % update fitting function
%     fitfun = 'funtaylor1';
%     PST.fitfun = 'funtaylor1';   
% end

%error('Stopping here to debug on 20150609\n');


% Calculate costs of null and initial models
switch ianneal
    case {0,1,2,3,4,5,6}
        % WARNING: this makes an assumption about the objective function
        % costs00 = rarcm(xd,zeros(size(xd)));
        % cost00 = nanmean(costs00)/DNPC;
        % evaluate costs of null model and initial model in cycles
        p00 = zeros(size(p0)); % null model
        
%         %cost00  = feval(objfun,p00,fitfun,DST,PST,TST); % cost of null model
%         fprintf(1,'Starting to calculate cost of null    estimate\n');
%         cost00  = feval(objfun,p00,'funnull',DST,PST00,TST); % cost of null model
%         fprintf(1,'Starting to calculate cost of initial estimate\n');
%         cost0   = feval(objfun,p0, fitfun,DST,PST0,TST); % cost of initial model
%         msig = nan(size(p0)); % uncertainty of model parameters
        % 20160524 reduce number of arguments to 4
       %cost00  = feval(objfun,p00,fitfun,DST,PST,TST); % cost of null model
        fprintf(1,'Starting to calculate cost of null    estimate\n');
        PST00.fitfun = 'funnull';
        objfun
        %cost00  = feval(objfun,DST,PST00,TST); % cost of null model
        cost00  = feval(str2func(objfun),DST,PST00,TST); % cost of null model
        fprintf(1,'Starting to calculate cost of initial estimate\n');
        PST0.fitfun = fitfun;
        %cost0   = feval(objfun,DST,PST0,TST); % cost of initial model
        cost0   = feval(str2func(objfun),DST,PST0,TST); % cost of initial model
        msig = nan(size(PST0.p0)); % uncertainty of model parameters
        istatcode = 0;
    case 7
        fprintf(1,'Starting to calculate cost of null    estimate\n');
        PST00.fitfun = 'funnull';
        cost00  = feval(str2func(objfun),PST00.p1,DST,PST00,TST); % cost of null model
        fprintf(1,'Starting to calculate cost of initial estimate\n');
        PST0.fitfun = fitfun;
        cost0   = feval(str2func(objfun),PST0.p1,DST,PST0,TST); % cost of initial model
        msig = nan(size(PST0.p0)); % uncertainty of model parameters
        istatcode = 0;
%     case 3
%         fprintf(1,'ianneal is 3, so runnning untested code in %s...\n',mfilename);
%         clear options;
%         options(1)=2;
%         options(2)=1;
%         %[p0,F,model,energy,count]=simann1(objfun,bounds,options,fitfun,DST,PST,TST)
%         [DST,PSTtemp,TST] = simann1(objfun,bounds,options,fitfun,DST,PST,TST)
%         mdl0 = DST.phamod;
% 
%         % cost in radians
%         %   2010JAN27 need to wrap model to prevent negative values of cost
%         cost00 = mean(rarcm(DST.phaobs,zeros(size(DST.phaobs))))/DNPC;
%         cost0  = mean(rarcm(DST.phaobs,           rwrapm(mdl0)))/DNPC;
%      case 4
%         fprintf(1,'ianneal is 4, so runnning untested code in %s...\n',mfilename);        
%         p00 = zeros(size(p0)); % null model
% %         cost00  = feval(objfun,p00,fitfun,DST,PST,TST); % cost of null model
% %         cost0   = feval(objfun,p0, fitfun,DST,PST,TST); % cost of initial model        
%         cost00  = feval(objfun,p00,fitfun,DST,PST,TSTP); % cost of null model
%         cost0   = feval(objfun,p0, fitfun,DST,PST,TSTP); % cost of initial model        
       
    otherwise
        error(sprintf('Unknown value of ianneal = %d\n',ianneal));
end


fprintf(1,'\n');
fprintf(1,'Total Average Cost of null     model = %12.4G %s for %10ld observations in data set for inversion\n',cost00,objlabel,ndata);
fprintf(1,'Total Average Cost of initial  model = %12.4G %s for %10ld observations in data set for inversion\n',cost0 ,objlabel,ndata);


switch ianneal
    case 0
        options = zeros(1,1);
    case {1, 2, 4, 5}
        % % options for simulated annealing, see help anneal
        %
        %      OPTIONS(1) = scale of cooling schedule (default = 4).  Higher numbers
        %                    produce more exhaustive searches.
        if pselect == 7 || ianneal == 4
            options(1) = 2;
        else
            options(1) = 4;
        end
        %options(1) = 20;
        %       OPTIONS(2) = number of individual annealing runs (default = 3).  Higher
        %                    numbers produce more exhaustive searches and reduce dependency
        %                    on correctly guessing critical temperature.
        if ianneal == 4
            options(2) = 1;  % only one run through SA for each iteration
        else
            options(2) = nsaruns;
        end
        %options(2) = 5;
        %       OPTIONS(3) = grid spacing (default = 4).  Higher numbers permit finer levels
        %                    of parameter discretization.
        if pselect == 7 || ianneal == 4
            options(3) = 2;
        else
            %options(3) = 10; % necessary for tests with simulated data
            options(3) = 4; % necessary for tests with simulated data
        end
        %       OPTIONS(4) = temperature scale. To tweak this parameter, inspect
        %                    a graph of 'energy'.  Values higher than 3 or lower than 1 are
        %                    seldom, if ever, warranted.  The default (0) tries several different
        %                    values and works well for most problems.
        %                    Smaller values lead to higher temperatures
        options(4) = 0;
        %options(4) = 1e-1;
        %       OPTIONS(5) = flag that tells the algorithm whether the objective function
        %                    can accept matrix input (default = 0); set to '1' if yes.  The
        %                    objective function should accept models stored columnwise and
        %                    return a vector of costs.  Writing the objective function this
        %                    way can increase algorithm speed by 15-20%.
        options(5) = 0;
        %     OPTIONS(6) = flag that tells the algorithm whether to try improve the solution:
        %                  0 = not try to improve (default)
        %                  1 = try the Nelder-Mead simplex method using
        %                      'fminsearch' function in MATLAB optimization toolbox
        %                  2 = try the Constrained optimization method using
        %                      'fmincon' function with in MATLAB optimization toolbox
        %                       and 'interior-point' algorithm
        %                  3 = evaluate paramater uncertainties using
        %                      jackknife analysis and fmincon option above
        options(6) = 0;
        %options(6) = saopt6;
        %       OPTIONS(7) = flag the tells the algorithm whether to display informative
        %                    output (default = 0); set to '1' if yes.
        %                    2 logs misfit and parameters values to file
        %options(7) = 1; % This option is set below
        %     OPTIONS(8) = flag to use MATLAB's distributed computing routines
        %                  set to the number of available processors [default is 1]
        options(8) = nprocessors;
        %     OPTIONS(9) = flag to Initialize random number generator.
        %                  0 == do not initialize [default]
        %                  1 == initialize once
        options(9) = 1;
    case 3
        clear options;
        options(1) = 1; % inverse problem
        options(2) = 0; % use files named PST.IN and DST.DAT
    case 6
        % Markov Chain Monte Carlo
        nskip = 1000;
        %niter = 1e8; % exhaustive
        niter = 1e6; % quick test
        nburnin = 1000;
        mstepscale = 0.1;
    case 7 
        % constrained optimization
        clear options;
    otherwise
        error(sprintf('Unknown value of ianneal = %d\n',ianneal));
end

switch ianneal
    case 0
        disp 'Skipping Annealing'
        p1 = p0;acosts1(1)=NaN;
        cost1 = cost0;
    case 1
        fprintf (1,'\nStarting Constrained Optimization...\n');
        options(7) = 1;acosts1(1)=NaN;
        tanneal=tic;
        %[p1,f,model,energy,count] = anneal5(objfun,bounds,options,fitfun,DST,PST,TST);
        
        %function [mhat,fval,model,energy,count]=constrainedopt1(FUN,bounds,OPTIONS,varargin)
        [p1,f,model,energy,count] = constrainedopt1(objfun,bounds,options,fitfun,DST,PST0,TST);
        fprintf (1,'\nAnnealing ended after %15.0f seconds\n',toc(tanneal));
        msig = nan(size(p0));      
    case 2
        %start annealing
        if idatatype1 == -1
            warning('Annealing on range gradient observable is not yet implemented.\n');
        end
        fprintf (1,'\nStarting annealing with recording....\n');
        options(7) = 2;  % use this for anneal2
        tanneal=tic;
        %[p1,f,trials,energy,count] = anneal5(objfun,bounds,options,fitfun,DST,PST0,TST);
        [p1,f,trials,energy,count] = anneal6(objfun,options,DST,PST0,TST);
        
        fprintf (1,'\nAnnealing ended after %15.0f seconds\n',toc(tanneal));
        trials = trials';
        acosts1= trials(:,1);        % cost values are in first column
        temps  = trials(:,2);        % temperatures are in second column
        trials = trials(:,3:end);% trial parameter values are in remaining columns
%         % calculate uncertainties in estimated parameters by bootstrap
%         % resampling (100 times slower than simulated annealing)
%         if saopt6 == 3
%             options(7) = 0;  % quiet
%            %options(7) = 2;  % verbose     
%            %[pj1,msig,costj1]=fminjackknife2(objfun,bounds,p1,options,fitfun,DST,PST,TST);
%             [pj1,msig,costj1,pboot]=bootstrapanneal1(objfun,bounds,p1,options,fitfun,DST,PST0,TST);
%             nf=nf+1;
%             h(nf) = plot_pairwise_correlations1(pboot,pnames);
%             feval(printfun,sprintf('%s_BOOTSCATTER',runname));
%         end
    case 3
        fprintf (1,'\nStarting annealing using SIMANN...\n');
        acosts1(1)=NaN;
        tanneal=tic;
        %[p1,F,model,energy,count]=simann1(objfun,bounds,options,fitfun,DST,PST,TST);
        [DST,PST1,TST] = simann1(objfun,bounds,options,fitfun,DST,PST0,TST)
        fprintf (1,'\nAnnealing using SIMANN ended after %15.0f seconds\n',toc(tanneal));
    case 4
        if surrogate ~= 1
           error(sprintf('surrogate (%d) is not set',surrogate));
        end
        imode = 2;
        
        maxiter = nsaruns;
       
        [TSTP, PST4, trials] = anneal_with_surrogate(objfun, DST, PST0, TST, maxiter...
            , verbose, imode, nf, runname, printfun, options);
        p1 = PST4.p1;
        msig = nan(size(p1));
        trials = trials';
        acosts1= trials(:,1);        % cost values are in first column
        temps  = trials(:,2);        % temperatures are in second column
        trials = trials(:,3:end);% trial parameter values are in remaining columns
    case 5
        fprintf (1,'\nStarting optimization using gridsearch...\n');
        acosts1(1)=NaN;
        options(7) = 2;  % verbose     
        tanneal=tic;
        %function [mhat,fval,model,energy,count,msig] = constrainedopt3(FUN,bounds,OPTIONS,fitfun,DST,PST,TST)
        %[p1,f,trials,energy,count,msig] = constrainedopt5(objfun,bounds,options,fitfun,DST,PST,TST);
        [p1,f,trials,energy,count,msig] = gridsearch1(objfun,bounds,options,fitfun,DST,PST,TST);
        fprintf (1,'\nOptimization ended after %15.2f seconds and %d evaluations\n',toc(tanneal),count);
    case 6                
        %[PSTout,h,mskip] = mcmc2(PST,DST,TST,mstepscale,nburnin,nskip,niter);
        [PST1,h,mskip] = mcmc3(PST0,DST,TST,mstepscale,nburnin,nskip,niter);
        p1 = PST1.p1;
        msig = PST1.sigma;
        acosts1(1)=NaN;
    case 7 
        % Constrained optimization
        %% Set constraints define matrix for negative values
%         Alessthan = eye(numel(p0));
%         blessthan = zeros(numel(p1),1);
        % Alessthan = [0, 0, -1, 0];
        % blessthan = 0;
        Alessthan = [];
        blessthan = [];
            
        
        %% Run constrained inversion
        %    X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
        %     the default optimization parameters replaced by values in OPTIONS, an
        %     argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
        %     details. For a list of options accepted by fmincon refer to the
        %     documentation.
        % options for algorithm
        options = optimoptions('fmincon','Display','Iter'); % run interior-point algorithm
        %funobj = 'funcoststdnres4'; % name of objective function
        objfun
        fh = str2func(objfun); % function handle
        p1 = fmincon(@(p1,pst,dst,tst) (fh(p1,DST,PST0,TST)),p0,Alessthan,blessthan,[],[],lb,ub,[],options);
        msig = PST0.sigma;
        acosts1(1)=NaN;
    otherwise
        error(sprintf('Unknown value of ianneal = %d\n',ianneal));
end; % switch on ianneal

%% save final estimate in structure PST1
PST1       = PST0;
PST1.p0    = colvec(p0); 
PST1.p1    = colvec(p1);
PST1.sigma = colvec(msig);



%% cost of final estimate
if ianneal == 7
    cost1 = feval(objfun,PST1.p1,DST,PST1,TST);
else
    cost1 = feval(objfun,DST,PST1,TST);
end

%% print results
fprintf(1        ,'\n\nI          Parameter_Name               P0(initial)     P1(final)     Adj   PlusMinusBnd\n');
for i = 1:mparam
   fmtstr = getfmt(p1(i),pnames{i});
   fprintf(1,fmtstr,pflag{i},i,pnames{i},p0(i),p1(i),p1(i)-p0(i),NaN,NaN,(ub(i)+lb(i))/2.0);
end


for fd=[1 fidtxtout]
   fprintf(fd,'\n');
   fprintf(fd,'Total Average Cost of null    model = %12.4G %s for %10ld observations in data set for inversion\n',cost00, objlabel,ndata);
   fprintf(fd,'Total Average Cost of initial model = %12.4G %s for %10ld observations in data set for inversion\n',cost0,  objlabel,ndata);
   fprintf(fd,'Total Average Cost of final   model = %12.4G %s for %10ld observations in data set for inversion\n',cost1,  objlabel,ndata);

   % Final is WORSE than initial
   if cost1 > cost0
      fprintf(fd,'WARNING: Final estimate is WORSE than initial estimate\n')
   end
end


% % make a 3-D plot using final estimate
% nf=nf+1;
% h(nf) = plot_obsmod3d(DST);
% feval(printfun,sprintf('%s_OBSMOD3D',runname));
% 
% % make plot of observed versus modeled
% nf=nf+1;h(nf)=figure;subplot(2,1,1);
% plot(rwrapm(DST.phaobs)/2.0/pi,'ro');hold on;
% plot(rwrapm(DST.phamod)/2.0/pi,'k-');
% title('final estimate');
% ylabel('phase (cycles)');
% legend('observed','modeled');
% subplot(2,1,2);
% plot(rwrapm(DST.phaobs-DST.phamod)/2.0/pi,'bo-');
% xlabel('pixel index');ylabel('phase (cycle)');legend('residual');
% feval(printfun,sprintf('%s_Taylor1Check',runname));

clear h;
save('gipht.mat');

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));
return


