% GIPHT_step3
% Perform statistics
% For stand alone executable version, replace feval with feval2
%
% 2009-APR-13 slight modifications to case with anneal == 2
% 2009-DEC-07
% 2010-MAR-22 change totcost[00,0,1] to cost[00,0,1]
%to save calls to fitfundostat = 1;
% 2011-JUL-25
% 2012-JUL-10 helene

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

clearvars
load
fidtxtout = fopen(txtoutname,'a');

En00 = NaN;
En0  = NaN;
En1  = NaN;
Pmu  = NaN;
P00  = NaN;
P0   = NaN;
Rbar00 = NaN;
Rbar0  = NaN;
Rbar1  = NaN;
Kappa00 = NaN;
Kappa0  = NaN;
Kappa1  = NaN;
crit69 = NaN;
cstddev1 = NaN;

% decide how to handle statistics
switch pselect
    case {1,2,3,5} % observable is phase
        % Assume phase residuals are distributed as Von Mises to calculate
        % critical value of cost at 69 percent confidence
        istatcode = 1;
    case 7 % observable is gradient
        if strcmp(objfun,'funcostrarcscaled') == 1
            % Assume residuals of gradient values are distributed as Gamma to calculate
            % critical value of cost at 69 percent confidence
            istatcode = 6;
        else
            % Assume deviations of gradient values are distributed as Generalized Pareto to calculate
            % critical value of cost at 69 percent confidence
            istatcode = 2;
        end
    otherwise
        % skip statistics
        istatcode = 0;
end

switch ianneal
    case 0
        % no inversion has been performed, so skip statistical analysis
        istatcode = 0;
    case {1,2}
        if saopt6 == 3
            istatcode = 5; % use Bootstrap instead of critical value
        end
        fprintf(1,'IANNEAL = %d so ISTATCODE changed to %d\n',ianneal,istatcode);
    case 6
        % use uncertainties from MCMC
        istatcode = 7;
    otherwise
        istatcode = 0;
        fprintf(1,'IANNEAL = %d so ISTATCODE defaults to %d\n',ianneal,istatcode);
end



% make a set of parameter estimates with all geophysical parameters set to zero
p00 = zeros(size(p0));

% make a set of parameter estimates with all the geophysical parameters set to zero
%[pt, px, py, pz, pb, pd, pg, pnames] = chopparamvector(p1,mparams,me);
%p2 = [pt; px; py; pz; pb; pd; zeros(numel(pg),1)];
% pindex=get_parameter_index('epoch',mparams,me,fitfun);pindex=setdiff(1:me,pindex);
% p2=p1;p2(pindex)=0.0;

% make a set of parameter estimates with all the nuisance parameters set to zero
% keep only pre-fit geophysical parameters
%[pt, px, py, pz, pb, pd, pg, pnames] = chopparamvector(p0,mparams,me);
%p3 = [pt; zeros(numel(px),1); zeros(numel(py),1); zeros(numel(pz),1); zeros(numel(pb),1); zeros(numel(pd),1); pg];

% make a set of parameter estimates with all the nuisance parameters set to zero
% keep only post-fit geophysical parameters
%[pt, px, py, pz, pb, pd, pg, pnames] = chopparamvector(p1,mparams,me);
%p4 = [pt; zeros(numel(px),1); zeros(numel(py),1); zeros(numel(pz),1); zeros(numel(pb),1); zeros(numel(pd),1); pg];

%%evaluate fitting function at current value of parameters
switch ianneal
    case {0,1,2,3,4,5,6}
        %% range change in radians
        % null model
        %mdl00 = feval(fitfun,DST,PST00,TST); % null
        mdl00 = feval('funnull',DST,PST00,TST); % null
        %mdl00 = mdl00 + mean_direction(DST.phaobs);
        % initial model
        mdl0  = feval(fitfun,DST,PST  ,TST); % initl
        % final model
        mdl1  = feval(fitfun,DST,PST1 ,TST); % final
        %     case 3
        %         fprintf(1,'ianneal is 3, so runnning untested code in %s...\n',mfilename);
        %         clear options;
        %         options(1)=2;
        %         options(2)=0;
        %         PST00.p0(49)
        %         PST.p0(49)
        %         PST1.p0(49)
        %         [DST,PST,TST]=simann1(objfun,bounds,options,fitfun,DST,PST00,TST);
        %         mdl00 = DST.phamod;
        %         disp 'mdl00'; sum(mdl00)
        %
        %         [DST,PST,TST]=simann1(objfun,bounds,options,fitfun,DST,PST,TST);
        %         mdl0 = DST.phamod;
        %         disp 'mdl0 '; sum(mdl0)
        %
        %         [DST,PST,TST]=simann1(objfun,bounds,options,fitfun,DST,PST1,TST);
        %         mdl1 = DST.phamod;
        %         disp 'mdl1'; sum(mdl1)
    otherwise
        error(sprintf('Unknown value of ianneal %d\n',ianneal));
end

% 2011-NOV-30 take real part
if isreal(mdl00) ~= 1
    warning(sprintf('Found complex values in mdl00. Max abs(imaginary) is %g\n',max(imag(mdl00))));
    mdl00 = real(mdl00);
end
if isreal(mdl0) ~= 1
    warning(sprintf('Found complex values in mdl0. Max abs(imaginary) is %g\n',max(imag(mdl0))));
    mdl0 = real(mdl0);
end
if isreal(mdl1) ~= 1
    warning(sprintf('Found complex values in mdl00. Max abs(imaginary) is %g\n',max(imag(mdl00))));
    mdl1 = real(mdl1);
end
% wrapped values
wrm00= rwrapm(mdl00); % radians [-pi, +pi]
wrm0 = rwrapm(mdl0);  % radians [-pi, +pi]
wrm1 = rwrapm(mdl1);  % radians [-pi, +pi]

% 20120507
% Wrap residuals for all observables because needed
% for phase plots in gipht_step5
% Unwrapped residuals are calculated in gipht_step4
res00=rwrapm(DST.phaobs-mdl00);  % radians
res0 =rwrapm(DST.phaobs-mdl0);   % radians
res1 =rwrapm(DST.phaobs-mdl1);   % radians

% Mean direction of residuals, allowing arbitrary mean direction
% mnd00= mean_direction(reshape(double(res00)/DNPC,numel(res00),1)*2.0*pi)/2.0/pi;  % in cycles
% mnd0 = mean_direction(reshape(double(res0) /DNPC,numel(res0) ,1)*2.0*pi)/2.0/pi;  % in cycles
% mnd1 = mean_direction(reshape(double(res1) /DNPC,numel(res1) ,1)*2.0*pi)/2.0/pi;  % in cycles
mnd00= mean_direction(res00*DNPC)/2.0/pi;  % in cycles
mnd0 = mean_direction( res0*DNPC)/2.0/pi;  % in cycles
mnd1 = mean_direction( res1*DNPC)/2.0/pi;  % in cycles

% % Circular Mean Deviations of residuals, allowing arbitrary mean direction
% cmd0 = circular_mean_deviation(reshape(double(res0)/DNPC,numel(res0),1)*2*pi,1)/2/pi;  % in cycles
% cmd1 = circular_mean_deviation(reshape(double(res1)/DNPC,numel(res1),1)*2*pi,1)/2/pi;  % in cycles
%
% % Circular Mean Deviations of residuals, setting mean direction to zero
% cmd0_0 = circular_mean_deviation(reshape(double(res0)/DNPC,numel(res0),1)*2*pi,0)/2/pi;  % in cycles
% cmd0_1 = circular_mean_deviation(reshape(double(res1)/DNPC,numel(res1),1)*2*pi,0)/2/pi;  % in cycles

% unwrapped versions
umd0 = mdl0 + double(res0);
umd1 = mdl1 + double(res1);

if strcmp(objfun,'funcostrarc') == 1
    costs00 = rarcm(DST.phaobs,wrm00);
    costs0  = rarcm(DST.phaobs,wrm0);
    costs1  = rarcm(DST.phaobs,wrm1);
elseif strcmp(objfun,'funcostrarcscaled') == 1
    costs00 = rarcm(DST.phaobs,wrm00) ./ DST.phasig;
    costs0  = rarcm(DST.phaobs,wrm0)  ./ DST.phasig;
    costs1  = rarcm(DST.phaobs,wrm1)  ./ DST.phasig;
else
    error(sprintf('Unknown value of objfun %s\n',objfun));
end
% values of total costs in cycles
cost00 = feval(objfun,p00,fitfun,DST,PST,TST);
cost0  = feval(objfun,p0, fitfun,DST,PST,TST);
cost1  = feval(objfun,p1, fitfun,DST,PST,TST);

% mean resultant lengths of resdiduals
Rbar00 = rbarrad(res00);
Rbar0  = rbarrad(res0);
Rbar1  = rbarrad(res1);

if cost1 > cost0
    fprintf(1,'WARNING: final estimate (%.4f cycle/datum) is NO BETTER than initial estimate (%.4f cycle/datmum).\n',cost1,cost0);
end
if isfinite(cost0) ~= 1 || isfinite(cost1) ~= 1
    fprintf(1,'WARNING: undefined cost for final estimate (%.4f cycle/datum) or initial estimate (%.4f cycle/datmum).\n',cost1,cost0);
    fprintf(1,'WARNING: not calculating statistics\n');
    istatcode = 0;
end

if istatcode ~= 0
    printstats= 1;
    %Circular standard deviation of final residuals
    %cstddev1  = circular_stddev(2*pi*reshape(res1,numel(res1),1))/2/pi;
    cstddev1  = circular_stddev(2*pi*reshape(double(res1)/DNPC,numel(res1),1))/2/pi;
    
    %     % make raw histogram of residuals
    %     nbins = floor(ndata/100);
    %     nbins = nbins - mod(nbins,2); % want an even number of bins
    %     if nbins < 10
    %         nbins = 10;
    %     end
    %     if nbins > 50
    %         nbins = 50;
    %     end
    %     nf=nf+1;h(nf)=figure;
    %     subplot(2,1,1);
    %     hist(colvec(res0)/DNPC,nbins);
    %     title('histogram of initial residuals');xlabel('residual (cycle)');ylabel('N occurrences');
    %     subplot(2,1,2);
    %     hist(colvec(res1)/DNPC,nbins);
    %     title('histogram of final residuals');xlabel('residual (cycle)');ylabel('N occurrences');
    %     feval(printfun,sprintf('%s_HIST2',runname));
    
    
    switch istatcode
        case 1  % Assume Von Mises distribution
            fprintf(1,'Evaluating Von Mises distribution for residual phase values.\n')
            
            % make QQ plot of residuals
            nf=nf+1;h(nf)=figure;
            qqplotvonmises(reshape(double(res1)/DNPC,numel(res1),1));
            feval(printfun,sprintf('%s_QQPLOT_von_mises',runname));
            
            % Test for von Miseness
            [Sm,VMnessString1] = vonmisesness (colvec(res1));    % Mardia and Jupp
            [U2,VMnessString2] = vonmisesness2(colvec(res1));    % Fisher
            
            if numel(strfind(VMnessString1,'ARE')) == 0 || numel(strfind(VMnessString2,'ARE')) > 0
                fprintf(1,'Assumption of von Mises Distribution is valid. ISTATCODE = %d\n',istatcode);
                %             % test if mean direction of residuals is zero
                %             [En00, Pmu] = testmu0(colvec(res00),0,0,0.05);
                %             [En0 , Pmu] = testmu0(colvec(res0) ,0,0,0.05);
                %             [En1,  Pmu] = testmu0(colvec(res1) ,0,0,0.05);
                %
                % test goodness of fit
                P00 = testkappas3(Rbar00,numel(DST.phaobs),Rbar1,numel(DST.phaobs));
                P0  = testkappas3(Rbar0, numel(DST.phaobs),Rbar1,numel(DST.phaobs));
                
                Rbar00 = rbarrad(res00);
                Rbar0  = rbarrad(res0);
                Rbar1  = rbarrad(res1);
                Kappa00 = batschelet_inv(Rbar00);
                Kappa0  = batschelet_inv(Rbar0);
                Kappa1  = batschelet_inv(Rbar1);
                
                if P0 < 0.025  % 1-sided
                    fprintf(1,'Final Model is significantly different (BETTER) than Initial Model with 95 percent confidence\n');
                elseif P0 > 0.975
                    fprintf(1,'Final Model is significantly different (worse!) than Initial Model with 95 percent confidence\n');
                else
                    fprintf(1,'Final Model is NOT significantly different than Initial Model with 95 percent confidence\n');
                end
                
                % Find critical levels of cost. The test statistic is normally distributed
                %Rbar68=testkappas_inv1(-1.00,ndata,Rbar1,ndata) % Z = 1.0  1-sigma = 69 percent confidence
                %Rbar95=testkappas_inv1(-1.96,ndata,Rbar1,ndata) % Z = 1.96           95
                % Find critical levels of cost. The test statistic is normally distributed
                %disp('Using testkappas_inv1')
                %Rbar68=testkappas_inv1(-0.48,ndata,Rbar1,ndata) % Z = -0.48  is P = 31% one sided
                %     Rbar95=testkappas_inv1(-1.65,ndata,Rbar1,ndata) % Z = -1.645 is P = 5 % one sided
                % % Find critical levels of cost. The test statistic is normally distributed
                %         disp('Usingtestkappas_inv2')
                %         Rbar68=testkappas_inv2( 0.31,ndata,Rbar1,ndata) % Z = -0.48  is P = 31% one sided
                %         Rbar95=testkappas_inv2( 0.05,ndata,Rbar1,ndata) % Z = -1.645 is P = 5 % one sided
                fprintf(1,'Calling testkappas_inv3 with Rbar1 = %.6f ndata = %d\n',Rbar1,ndata)
                Rbar69 = testkappas_inv3(0.31,ndata,Rbar1,ndata) % Z = -0.48  is P = 31% one sided
                
                %     if Rbar69 > Rbar1
                %         warning(sprintf('Rbar69 (%.4f) is greater than Rbar1 (%.4f)\n',Rbar69,Rbar1));
                %     end
                %     %     if Rbar95 > Rbar1
                %     %         warning(sprintf('Rbar95 (%.4f) is greater than Rbar1 (%.4f)\n',Rbar95,Rbar1));
                %     %     end
                Kappa69=batschelet_inv(Rbar69) % critical value of concentration parameter
                %   Kappa95=batschelet_inv(Rbar95);
                %    crit69=kappa2cmd(Kappa69,ndata,mean_direction(reshape(res1,numel(res1),1))*pi*2)/pi/2;  % critical value of circular mean deviation (cost) in cycles
                %    crit95=kappa2cmd(Kappa95,ndata,mean_direction(reshape(res1,numel(res1),1))*pi*2)/pi/2;  % critical value of circular mean deviation (cost) in cycles
                
                %crit69=kappa2cmd(Kappa69,ndata,mean_direction(reshape(double(res1)/DNPC,numel(res1),1))*pi*2)/pi/2;  % critical value of circular mean deviation (cost) in cycles
                %crit69=kappa2cmd(Kappa69,ndata,mean_direction(colvec(res1)))/2.0/pi  % critical value of circular mean deviation (cost) in cycles
                crit69=kappa2cmd(Kappa69,ndata)/2.0/pi  % critical value of circular mean deviation (cost) in cycles, assuming zero mean
                %crit95=kappa2cmd(Kappa95,ndata,mean_direction(reshape(double(res1)/DNPC,numel(res1),1))*pi*2)/pi/2;  % critical value of circular mean deviation (cost) in cycles
                %fprintf(1,'Critical Value of Circ. Mean Dev. = %.4f cycles per datum for %6d observations at 95 percent confidence\n',crit95, numel(xd));
                %crit69=kappa2cost(Kappa69);  % critical value of cost in cycles
                %crit95=kappa2cost(Kappa95);  % critical value of cost in cycles
                %fprintf(1,'Critical Value of cost = %.4f cycles per datum for %6d observations at 69 percent confidence\n',crit69, numel(xd));
                %fprintf(1,'Critical Value of cost = %.4f cycles per datum for %6d observations at 95 percent confidence\n',crit95, numel(xd));
            else
                istatcode = -1 * istatcode;
            end
            
            
        case 2
            fprintf(1,'Evaluating Generalized Pareto distribution for angular deviation between observed and modeled values.\n')
            
            % estimate parameters in distribution
            dist = 'Generalized Pareto';
            %
            %             % critical value
            crit69 = test_generalized_paretos(costs1/DNPC,cost1)
            
            fprintf(1,'Evaluating %s distribution for angular deviations between observed and modeled.\n',dist)
            [h,phat,chi2gof_h,chi2gof_p,chi2gof_stats] = qqplot(colvec(costs1)/DNPC,dist);
            
            feval(printfun,sprintf('%s_QQPLOTdeviations_%s'   ,runname,strrep(dist,' ','_')),h(1));
            feval(printfun,sprintf('%s_HISTOGRAMdeviations_%s',runname,strrep(dist,' ','_')),h(2));
            
            if chi2gof_h == 0
                istatcode = abs(istatcode);
                fprintf(1,'Assumption of %s Distribution is valid. ISTATCODE = %d\n',dist,istatcode);
            else
                istatcode = -1*abs(istatcode);
                fprintf(1,'Assumption of %s Distribution is NOT valid. ISTATCODE = %d\n',dist,istatcode);
            end
        case 6
            dist = 'Gamma';
            fprintf(1,'Evaluating %s distribution for scaled angular deviations between observed and modeled.\n',dist)
            [h,phat,chi2gof_h,chi2gof_p,chi2gof_stats] = qqplot(colvec(costs1)/DNPC,dist);
            feval(printfun,sprintf('%s_QQPLOTdeviations_%s'   ,runname,strrep(dist,' ','_')),h(1));
            feval(printfun,sprintf('%s_HISTOGRAMdeviations_%s',runname,strrep(dist,' ','_')),h(2));
            
            if chi2gof_h == 0
                istatcode = abs(istatcode);
                fprintf(1,'Assumption of %s Distribution is valid. ISTATCODE = %d\n',dist,istatcode);
                crit69 = test_gammas(colvec(costs1)/DNPC);
            else
                istatcode = -1*abs(istatcode);
                fprintf(1,'Assumption of %s Distribution is NOT valid. ISTATCODE = %d\n',dist,istatcode);
                crit69 = NaN;
            end
        case 5
            crit69 = NaN;
            fprintf(1,'Using Bootstrap statistics ISTATCODE = %d\n',istatcode);
        case 7
            crit69 = NaN;
            fprintf(1,'Using MCMC statistics ISTATCODE = %d\n',istatcode);
        otherwise
            warning(sprintf('Unknown value of ISTATCODE = %d\n',istatcode));
    end
    
    if ismember(abs(istatcode),[0,1,2,5,6]) == 1
        psig = nan(size(p0));
        if isfinite(crit69) == 0
            warning(sprintf('Critical value is not defined.\nApproximating critical value using non-parametric quantiles.\n'));
            if isfinite(acosts1(1)) == 1
                crit69 = quantile(acosts1(isfinite(acosts1)),0.69);
                istatcode = 3;
                fprintf(1,'Approximating critical value using non-parametric quantiles . ISTATCODE = %d\n',istatcode);
            else
                %crit69a = quantile(costs1/DNPC,0.69);
                % 20130701
                iok=find(isfinite(costs1)==1);
                crit69a = quantile(costs1(iok)/DNPC,0.69);
                crit69b = cost1 + quantile(costs1/DNPC,0.69) - nanmin(colvec(costs1/DNPC));
                crit69 = nanmin([crit69, crit69a, crit69b]);
                istatcode = 4;
                fprintf(1,'Approximating critical value using non-parametric quantiles on residuals. ISTATCODE = %d\n',istatcode);
            end
        end
        if crit69 <= cost1 || isfinite(crit69) == 0
            warning(sprintf('Critical value is lower than final (minimum) value.\nApproximating critical value using empirical values.\n'));
            crit69 = cost1 + (1.0 - 0.884322217359915) * abs(cost1 - cost0);
            istatcode = 6;
            fprintf(1,'Approximating critical value using empirical values. ISTATCODE = %d\n',istatcode);
        end
        
        if istatcode ~= 0 && isfinite(crit69) == 1
            % evaluate cost for each parameter
            nslices = 100;
            %p1vals = zeros(numel(acosts1),mparam); % initialize  this matrix
            psig = nan(mparam,1);            % intialize sigmas
            lsig = nan(mparam,1);            % intialize sigmas
            rsig = nan(mparam,1);            % intialize sigmas
            ymin = cost1  - 0.05*abs(cost0-cost1);
            if isfinite(acosts1(1)) == 1
                %ymax = max(acosts1); % show costs of all trials
                ymax = quantile(acosts1,0.69); % show only the best
                ymax = nanmax([ymax cost1]);
                ymax = nanmax([ymax (crit69 + 0.05*abs(crit69-cost1))]);
                ymax = nanmax([ymax (cost0  + 0.05*abs(cost0 -cost1))]);
            else
                ymax = cost0  + 2.0*abs(cost0-cost1);
                %ymax = 1.1*cost0;
            end
            costdiffthreshl = 5*abs(cost0-cost1)/nslices;
            costdiffthreshr = 5*abs(cost0-cost1)/nslices;
            adjthresh = 1.0e-5;
            for i=1:mparam
                %for i = 58;
                lb=PST1.lb;
                ub=PST1.ub;
                p0=PST1.p0;
                p1=PST1.p1;
                lsig(i) = lb(i);leftisset=0;   % left  1-sigma limit
                rsig(i) = ub(i);rigtisset=0;   % right 1-sigma limit
                % remove underscore to avoid problem with axis label
                tmpstr = char(pnames{i});
                pname_no_=strrep(tmpstr,'_','\_');
                
                db = (ub(i)-lb(i));
                nbins = floor(1.1*nslices);
                costs5l = nan(nbins,1);
                costs5r = nan(nbins,1);
                costs5lp = nan(nbins,1);
                costs5rp = nan(nbins,1);
                dpl = nan(nbins,1);
                dpr = nan(nbins,1);
                if db > abs(p1(i)-p0(i)) && numel(strfind(pnames{i},'Offset')) == 0
                    fprintf(1,'Varying parameter: %s\n',pnames{i});
                    if isfinite(acosts1(1)) == 1
                        p1vals = trials(:,i)';
                    end
                    dp = p1(i) + db * [-2 2];
                    
                    % find left side of 69 percent confidence interval
                    if istatcode == 5 % use Bootstrap uncertainty
                        lsig(i) = p1(i) - msig(i);
                        leftisset = 1;
                    elseif ianneal == 2 && isfinite(crit69) == 1
                        % find minimal trial value, per Tabrez' suggestion
                        % 2012-JAN-03
                        isubcritl = find(acosts1 <= crit69);
                        if numel(isubcritl) > 0
                            lsig(i) = nanmin(trials(isubcritl,i));
                        else
                            warning('No cost values less than critical value.');
                            lsig(i) = lb(i);
                        end
                        leftisset = 1;
                    else
                        warning('Cannot calculate uncertainty for estimated parameter');
                        lsig(i) = NaN;
                    end
                    
                    % find right side of 69 percent confidence interval
                    if istatcode == 5 % use Bootstrap uncertainty
                        rsig(i) = p1(i) + msig(i);
                        rigtisset = 1;
                    elseif ianneal == 2 && isfinite(crit69) == 1
                        isubcritr = find(acosts1 <= crit69);
                        if numel(isubcritr) > 0
                            rsig(i) = nanmax(trials(isubcritr,i));
                        else
                            warning('No cost values less than critical value.');
                            rsig(i) = ub(i);
                        end
                        
                        rigtisset = 1;
                    else
                        warning('Cannot calculate uncertainty for estimated parameter');
                        rsig(i) = NaN;
                    end
                    
                    
                    
                    % find half-width of 69 percent confidence interval as average
                    psig(i) = (abs(lsig(i)-p1(i)) + abs(rsig(i)-p1(i)) )/2;   %
                    fprintf(1,'%32s %20.10e +/- %20.10e\n',char(pnames{i}),p1(i),psig(i));
                    
                    %if db > abs(p1(i)-p0(i)) && numel(strfind(pnames{i},'Offset')) == 0
                    nf=nf+1;h(nf)=figure;
                    axis([-Inf,+Inf,ymin,ymax]);
                    hold on;
                    
                    if isfinite(acosts1(1)) == 1
                        plot(p1vals,acosts1,'k.');
                    end
                    
                    %               dp = [dpl dpr];
                    %costs5 = [costs5l costs5r];
                    iok = find(isfinite(costs5l));
                    plot(dpl(iok),costs5l(iok),'Color','b','LineStyle','-','LineWidth',2);
                    iok = find(isfinite(costs5r));
                    plot(dpr(iok),costs5r(iok),'Color','b','LineStyle','-','LineWidth',2);
                    iok = intersect(find(isfinite(costs5lp))...
                        ,find(costs5lp > cost1));
                    %                 iok = intersect(iok,find(costs5lp < ymax));
                    %                 iok = intersect(iok,find(dpl>lb(i)));
                    plot(dpl(iok),costs5lp(iok),'Color','g','LineStyle','-','LineWidth',2);
                    iok = intersect(find(isfinite(costs5rp))...
                        ,find(costs5rp > cost1));
                    %                 iok = intersect(iok,find(costs5rp < ymax));
                    %                 iok = intersect(iok,find(dpr<ub(i)));
                    
                    plot(dpr(iok),costs5rp(iok),'Color','g','LineStyle','-','LineWidth',2);
                    
                    plot ([lb(i) lb(i)],    [ymin  ymax],   'k:' ,'Clip','off','LineWidth',2); % lower bound
                    plot ([ub(i) ub(i)],    [ymin  ymax],   'k:' ,'Clip','off','LineWidth',2); % upper bound
                    plot ([p1(i) p1(i)],    [cost1 ymax],  'r-' ,'Clip','off','LineWidth',2);
                    plot ([p1(i) p1(i)],    [cost1 cost1], 'ro' ,'Clip','off','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',10);
                    plot ([p0(i) p0(i)],    [cost0 ymax],  'k--','Clip','off','LineWidth',2);
                    plot ([p0(i) p0(i)],    [cost0 cost0], 'ks' ,'Clip','off','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerSize',10);
                    plot ([min(min(dp)) max(max(dp))],   [cost0  cost0],             'k--','Clip','off','LineWidth',2);  % cost of initial model
                    plot ([min(min(dp)) max(max(dp))],   [crit69 crit69],            'k--','Clip','off','LineWidth',2);
                    plot ([p1(i)-psig(i) p1(i)+psig(i)], [cost1  cost1],'r-' ,'Clip','off','LineWidth',4); % 1-sigma error bar at base
                    h2 = text(max(max(dp))-0.20*(max(max(dp))-min(min(dp))),crit69,sprintf('69%%'));     set(h2,'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0],'HorizontalAlignment','Right','FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    h2 = text(min(min(dp))+0.10*(max(max(dp))-min(min(dp))),cost0, sprintf('initial'));  set(h2,'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0],'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    h2 = text(min(min(dp))+0.20*(max(max(dp))-min(min(dp))),cost1, sprintf('final'));    set(h2,'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[1 0 0],'Color',[1 0 0],'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    
                    
                    
                    if cost1 > 0.01
                        fixlabels(pname_no_,'','circular mean deviation (cycles/datum)','%10.4f');
                    else
                        fixlabels(pname_no_,'','circular mean deviation (cycles/datum)','');
                    end
                    
                    %                 h2=xlabel(pname_no_);set(h2,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    %                 h2=ylabel('circular mean deviation (cycles/datum)');set(h2,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    %                 %title(sprintf('%s %12.4G +/- %12.4G',pname_no_,p1(i),psig(i)));
                    outfmt=getfmt(p1(i),pname_no_);outfmt=outfmt(1:findstr('#',outfmt)+18);
                    h2=title(sprintf(outfmt,'P#',i,pname_no_,p1(i),psig(i)));set(h2,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
                    feval(printfun,sprintf('%s_param%03d',runname,i));
                end
            end
        end
        if cost1 > crit69
            fprintf(1,'WARNING: Final value of cost (%.4f) is higher (WORSE) than 69-percent critical value (%.4f).\n',cost1,crit69);
            fprintf(1,'WARNING: Skipping statistics\n');
            printstats = 0;
        end
    end
end


%fclose(fidtxtout);
if fidtxtout < 100
    fidtxtout=fopen(txtoutname,'a');
end
for fd=[1 fidtxtout]
    fprintf(fd,'Cost  of null  model   = %.7f cycles per datum for %6d observations in inverted data set %s\n',cost00,  ndata, runname);
    fprintf(fd,'Cost  of initl model   = %.7f cycles per datum for %6d observations in inverted data set %s\n',cost0,   ndata, runname);
    fprintf(fd,'Cost  of final model   = %.7f cycles per datum for %6d observations in inverted data set %s\n',cost1,   ndata, runname);
    fprintf(fd,'Cost  improvement      = %.7f cycles per datum for %6d observations in inverted data set %s\n',cost0-cost1,   ndata, runname);
if isfinite(crit69)==1
    fprintf(fd,'Critical value of cost = %.7f cycles per datum for %6d observations in inverted data set %s\n',crit69,  ndata, runname);
end
    if istatcode == 1
        
        fprintf(fd,'Mean direction of residuals from null     model              = %.4f cycles\n',mnd00);
        fprintf(fd,'Mean direction of residuals from initial  model              = %.4f cycles\n',mnd0);
        fprintf(fd,'Mean direction of residuals from final    model              = %.4f cycles\n',mnd1);
        
        fprintf(fd,'Circular standard deviation of residuals from final    model = %.4f cycles\n',cstddev1);
        
        fprintf(fd,'Mean Resultant Length for residuals for NULL  model Rbar00   = %10.4f\n',Rbar00);
        fprintf(fd,'Mean Resultant Length for residuals for initl model Rbar0    = %10.4f\n',Rbar0);
        fprintf(fd,'Mean Resultant Length for residuals for final model Rbar1    = %10.4f\n',Rbar1);
        
        fprintf(fd,'Concentration Parameter of residuals for NULL  model Kappa00 = %10.4f\n',Kappa00);
        fprintf(fd,'Concentration Parameter of residuals for initl model Kappa0  = %10.4f\n',Kappa0);
        fprintf(fd,'Concentration Parameter of residuals for final model Kappa0  = %10.4f\n',Kappa1);
        
        %         fprintf(fd,'Test Statistic distribututed as N(0,1) eta00                = %10.4f\n',eta00);
        %         fprintf(fd,'Test Statistic distribututed as N(0,1) eta0                 = %10.4f\n',eta0);
        fprintf(fd,'P-value of Test Statistic P00                                = %10.4f\n',P00);
        fprintf(fd,'P-value of Test Statistic P0                                 = %10.4f\n',P0);
        
        fprintf(fd,'Test Statistic for mean direction of residuals En00          = %10.4f\n',En00);
        fprintf(fd,'Test Statistic for mean direction of residuals En0           = %10.4f\n',En0);
        fprintf(fd,'Test Statistic for mean direction of residuals En1           = %10.4f\n',En1);
    end
end

% 2012-JAN-11: this cannot work because i1 and i2 are not yet defined
% if printstats== 1
%     % print statistics broken down by pairs
%     titlestrs = cell(np);
%     for i=1:np
%         titlestrs{i} = sprintf('Pair %3d orbs %7d %7d years %6.1f to %6.1f Dt = %.4f yr Cost0 = %6.4f Cost1=%6.4f %s'...
%             ,i,iuniqorbs(kmast),iuniqorbs(kslav)...
%             ,tepochs(kmast),tepochs(kslav),tepochs(kslav)-tepochs(kmast)...
%             ,nanmean(double(costs0(i1:i2))/DNPC),nanmean(double(costs1(i1:i2))/DNPC) ...
%             ,runname);
%         for fid = [1 fidtxtout]
%             fprintf(fid,'%s %s\n',titlestrs{i},pfnames{i});
%         end
%     end
% end


% PRINT PARAMETERS TO SCREEN AND FILE
fprintf(1         ,'I          Name     P0(pre-fit)              P1(post-fit)     Adjust.   Sigma[%1d] Significance PlusMinusBound\n',istatcode);
fprintf(fidtxtout ,'I          Name     P0(pre-fit)              P1(post-fit)     Adjust.   Sigma[%1d] Significance PlusMinusBound\n',istatcode);
pflag = cell(mparam,1);
for i = 1:mparam
    adj = p1(i)-p0(i);
    if psig(i) > 0
        sadj = abs(adj/psig(i));
        pflag{i} = 'E#';
    else
        sadj = NaN;
        pflag{i} = 'F#';
    end
    
    outfmt = getfmt(p1(i),pnames{i});
    fprintf(1        ,outfmt,pflag{i},i,pnames{i} ,p0(i),p1(i),adj,psig(i),sadj,(ub(i)-lb(i))/2.0);
    fprintf(fidtxtout,outfmt,pflag{i},i,pnames{i}, p0(i),p1(i),adj,psig(i),sadj,(ub(i)-lb(i))/2.0);
end

% calculate some derived parameters
iq = 0;

% total volume change
iq = iq+1;
q0(iq) = (p0(get_parameter_index('Okada1_Length'  ,pnames))  ...
    *p0(get_parameter_index('Okada1_Width'  , pnames))  ...
    *p0(get_parameter_index('Okada1_Tensile', pnames))) ...
    + ( p0(get_parameter_index('Okada2_Length'  ,pnames))  ...
    *p0(get_parameter_index('Okada2_Width'  , pnames))  ...
    *p0(get_parameter_index('Okada2_Tensile', pnames))) ...
    +   p0(get_parameter_index('Mogi1_Volume_'  ,pnames))  ...
    +   p0(get_parameter_index('Mogi2_Volume_'  ,pnames));
q1(iq) = ( p1(get_parameter_index('Okada1_Length'  ,pnames))  ...
    *p1(get_parameter_index('Okada1_Width'  , pnames))  ...
    *p1(get_parameter_index('Okada1_Tensile', pnames))) ...
    + ( p1(get_parameter_index('Okada2_Length'  ,pnames))  ...
    *p1(get_parameter_index('Okada2_Width'  , pnames))  ...
    *p1(get_parameter_index('Okada2_Tensile', pnames))) ...
    +   p1(get_parameter_index('Mogi1_Volume_'  ,pnames))  ...
    +   p1(get_parameter_index('Mogi2_Volume_'  ,pnames));

% linearized propagation of uncertainties http://en.wikipedia.org/wiki/Propagation_of_uncertainty
% 2012-OCT-04
%qsig(iq) = 0;

tva = p1(get_parameter_index('Okada1_Length'  ,pnames)) ;
tvb = p1(get_parameter_index('Okada1_Width'  , pnames)) ;
tvc = p1(get_parameter_index('Okada1_Tensile', pnames)) ;
tsa = psig(get_parameter_index('Okada1_Length'  ,pnames)) ;
tsb = psig(get_parameter_index('Okada1_Width'  , pnames)) ;
tsc = psig(get_parameter_index('Okada1_Tensile', pnames)) ;

% correct the calculation of error below
if isfinite(tva*tvb*tvc*tsa*tsb*tsc) == 1
    qsig(iq) = (tvb*tvc*tsa)^2+(tva*tvc*tsb)^2+(tva*tvb*tsc)^2;
else
    % 2012-OCT-04
    qsig(iq) = NaN;
end


% Kurt 2012 JUL 10
if  isfinite(psig(get_parameter_index('YangPS_semimajor_axis_a__m______',pnames)))==1 ...
        && isfinite(psig(get_parameter_index('YangPS_semiminor_axis_b_in_m____',pnames)))==1
    qsig(iq) = NaN;
end
qsig(iq) = sqrt(qsig(iq));
qnames{iq} = sprintf('Total_Net_Volume_Increase_in_m3');


% DERIVED PARAMETERS FOR OKADA1
iii=get_parameter_index('Okada1_Length_in_m______________',pnames);
%if psig(iii) > 0
if p0(iii) > 0
    for i=1:10
        utmzone10(i,:)=utmzone0;
    end
%     % Geographic coordinates of Okada centroid
%     [Xcorners10,Ycorners10,Hcorners10,Ncorners] = disloc_to_seismo(p0(iii:iii+9));
%     [LatCorners10,LonCorners10]=utm2deg(Xcorners10,Ycorners10,utmzone10);
%     [Xcorners11,Ycorners11,Hcorners11,Ncorners] = disloc_to_seismo(p1(iii:iii+9));
%     [LatCorners11,LonCorners11]=utm2deg(Xcorners11,Ycorners11,utmzone10);
%     clear utmzone10;
%     
%     iq = iq+1;
%     q0(iq) = Xcorners10(10);
%     q1(iq) = Xcorners11(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada1_Centroid_Easting_in_m____');
%     iq = iq+1;
%     q0(iq) = Ycorners10(10);
%     q1(iq) = Ycorners11(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada1_Centroid_Northing_in_m___');
%     iq = iq+1;
%     q0(iq) = LatCorners10(10);
%     q1(iq) = LatCorners11(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada1_Centroid_latitude__in_deg');
%     iq = iq+1;
%     q0(iq) = LonCorners10(10);
%     q1(iq) = LonCorners11(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada1_Centroid_longitude_in_deg');
%     
%     %   fprintf(fd,'Initial UTM coordinates of 4 corners, upper center, and Centroid of Okada\n');
%     %   for i=[1 2 3 4 7 10]
%     %      fprintf(fd,'%s %12.4f %12.4f %12.4f\n',Ncorners{i},Xcorners0(i),Ycorners0(i),Hcorners0(i));
%     %   end
%     %   fprintf(fd,'Final   UTM coordinates of 4 corners, upper center, and Centroid of Okada\n');
%     %   for i=[1 2 3 4 7 10]
%     %      fprintf(fd,'%s %12.4f %12.4f %12.4f\n',Ncorners{i},Xcorners1(i),Ycorners1(i),Hcorners1(i));
%     %   end
%     %
%     %   fprintf(fd,'Initial Lon Lat of 4 corners, upper center, Centroid of Okada\n');
%     %   for i=[1 2 3 4 7 10]
%     %      fprintf(fd,'%s %12.4f %12.4f %12.4f\n',Ncorners{i},LonCorners0(i),LatCorners0(i),Hcorners0(i));
%     %   end
%     %   fprintf(fd,'Final   Lon Lat of 4 corners, upper center, Centroid of Okada\n');
%     %%   for i=[1 2 3 4 7 10]
%     %      fprintf(fd,'%s %12.4f %12.4f %12.4f\n',Ncorners{i},LonCorners1(i),LatCorners1(i),Hcorners1(i));
%     %   end
%     
%     
%     % conventional strike
%     iq=iq+1;
%     str = 180 + p0(get_parameter_index('Okada1_Strike',pnames));
%     if str > 360
%         str = str - 180;
%     end
%     q0(iq) = str;
%     str = 180 + p1(get_parameter_index('Okada1_Strike',pnames));
%     if str > 360
%         str = str - 180;
%     end
%     q1(iq) = str;
%     qsig(iq) = psig(get_parameter_index('Okada1_Strike',pnames));
%     qnames{iq} = sprintf('Okada1_Convt_strike_in_deg_CW_N');
%     
%     % Okada U1
%     % q0(9) = qg0(12);
%     % q1(9) = qg1(12);
%     %q0(9) = -1*p0(get_parameter_index('Okada1_RL_Strike_Slip',pnames));
%     %q1(9) = -1*p1(get_parameter_index('Okada1_RL_Strike_Slip',pnames));
%     %qnames{9} = sprintf('Okada_slip_U1_in_meters________');
%     % Okada U2
%     % q0(10) = -1*qg0(13);
%     % q1(10) = -1*qg1(13);
%     %q0(10) = -1*p0(get_parameter_index('Okada1_Downdip_Slip',pnames));
%     %q1(10) = -1*p1(get_parameter_index('Okada1_Downdip_Slip',pnames));
%     %qnames{10} =sprintf('Okada_slip_U2_in_meters________');
%     % Okada U3
%     % q0(11) = qg0(14);
%     % q1(11) = qg1(14);
%     %q0(11) = p0(get_parameter_index('Okada1_Tensile_Opening',pnames));
%     %q1(11) = p1(get_parameter_index('Okada1_Tensile_Opening',pnames));
%     %qnames{11} =sprintf('Okada_slip_U3_in_meters________');
    
%     %rake
%     iq = iq+1;
%     q0(iq) = 180*(atan2(-1*p0(get_parameter_index('Okada1_RL_Strike_Slip',pnames))...
%         ,p0(get_parameter_index('Okada1_Downdip_Slip'  ,pnames))...
%         ))/pi;  % change sign of U2 for negative dip bug in disloc
%     q1(iq) = 180*(atan2(-1*p1(get_parameter_index('Okada1_RL_Strike_Slip',pnames))...
%         ,p1(get_parameter_index('Okada1_Downdip_Slip',  pnames))...
%         ))/pi;  % change sign of U2 for negative dip bug in disloc
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Derived_Okada1_rake_in_deg_CCW_');
%     
%     % geometric potency
%     iq = iq+1;
%     q0(iq) = p0(get_parameter_index('Okada1_Length',pnames))...
%         *p0(get_parameter_index('Okada1_Width', pnames))...
%         *norm([p0(get_parameter_index('Okada1_RL_Strike_Slip',pnames))...
%         ,         p0(get_parameter_index('Okada1_Downdip_Slip',  pnames))...
%         ,         p0(get_parameter_index('Okada1_Tensile',       pnames))]);
%     q1(iq) = p1(get_parameter_index('Okada1_Length',pnames))...
%         *p1(get_parameter_index('Okada1_Width', pnames))...
%         *norm([p1(get_parameter_index('Okada1_RL_Strike_Slip',pnames))...
%         ,         p1(get_parameter_index('Okada1_Downdip_Slip'  ,pnames))...
%         ,         p1(get_parameter_index('Okada1_Tensile',       pnames))]);
%     qnames{iq} = sprintf('Derived_Okada1_potency_in_m3___');
    % qsig(iq) = NaN;
    %     qsig(iq) = 0;
    %     if isfinite(psig(get_parameter_index('Okada1_Tensile',pnames))) == 1
    %         qsig(iq) = qsig(iq) ...
    %             +((psig(get_parameter_index('Okada1_Length'  ,pnames))  ...
    %             * psig(get_parameter_index('Okada1_Width'  , pnames))   ...
    %             * psig(get_parameter_index('Okada1_Tensile', pnames))))^2;
    %   20130419 uncertainty is sqrt of sum of squares
%     qsig(iq) = sqrt((psig(get_parameter_index('Okada1_RL_Strike_Slip'  ,pnames)))^2  ...
%         + (psig(get_parameter_index('Okada1_Downdip_Slip'  , pnames)))^2   ...
%         + (psig(get_parameter_index('Okada1_Tensile', pnames)))^2);
%     2014-01-13 above gives strange results
    qsig(iq) = NaN;

end

% moment, assuming shear modulus = 30 GPa
iq = iq+1;
if q0(iq-1) > 0; q0(iq) = 3e10*q0(iq-1); else q0(iq) = NaN; end;
if q1(iq-1) > 0; q1(iq) = 3e10*q1(iq-1); else q1(iq) = NaN; end;
qnames{iq} = sprintf('Derived_Okada1_moment_in_Nm____');
qsig(iq) = NaN;

% magnitude
iq = iq+1;
if q0(iq-1) > 0; q0(iq) = 2*log10(q0(iq-1))/3 - 6.03; else q0(iq) = NaN; end
if q1(iq-1) > 0; q1(iq) = 2*log10(q1(iq-1))/3 - 6.03; else q1(iq) = NaN; end
qnames{iq} = sprintf('Derived_Okada1_Mw______________');
qsig(iq) = NaN;

%end

% DERIVED PARAMETERS FOR OKADA2
iii=get_parameter_index('Okada2_Length_in_m______________',pnames);
%if psig(iii) > 0
if p0(iii) > 0
    for i=1:10
        utmzone10(i,:)=utmzone0;
    end
    
    [Xcorners20,Ycorners20,Hcorners20,Ncorners] = disloc_to_seismo(p0(iii:iii+9));
    [LatCorners20,LonCorners20]=utm2deg(Xcorners20,Ycorners20,utmzone10);
    [Xcorners21,Ycorners21,Hcorners21,Ncorners] = disloc_to_seismo(p1(iii:iii+9));
    [LatCorners21,LonCorners21]=utm2deg(Xcorners21,Ycorners21,utmzone10);
    clear utmzone10;
%     
%     
%     iq = iq+1;
%     q0(iq) = Xcorners20(10);
%     q1(iq) = Xcorners21(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada2_Centroid_Easting_in_m____');
%     iq = iq+1;
%     q0(iq) = Ycorners20(10);
%     q1(iq) = Ycorners21(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada2_Centroid_Northing_in_m___');
%     iq = iq+1;
%     q0(iq) = LatCorners20(10);
%     q1(iq) = LatCorners21(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada2_Centroid_latitude__in_deg');
%     iq = iq+1;
%     q0(iq) = LonCorners20(10);
%     q1(iq) = LonCorners21(10);
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Okada2_Centroid_longitude_in_deg');
%     
%     % conventional strike
%     iq=iq+1;
%     str = 180 + p0(get_parameter_index('Okada2_Strike',pnames));
%     if str > 360
%         str = str - 180;
%     end
%     q0(iq) = str;
%     str = 180 + p1(get_parameter_index('Okada2_Strike',pnames));
%     if str > 360
%         str = str - 180;
%     end
%     q1(iq) = str;
%     qsig(iq) = psig(get_parameter_index('Okada2_Strike',pnames));
%     qnames{iq} = sprintf('Okada2_Convt_strike_in_deg_CW_N');
    
%     %rake
%     iq = iq+1;
%     q0(iq) = 180*(atan2(-1*p0(get_parameter_index('Okada2_RL_Strike_Slip',pnames))...
%         ,p0(get_parameter_index('Okada2_Downdip_Slip'  ,pnames))...
%         ))/pi;  % change sign of U2 for negative dip bug in disloc
%     q1(iq) = 180*(atan2(-1*p1(get_parameter_index('Okada2_RL_Strike_Slip',pnames))...
%         ,p1(get_parameter_index('Okada2_Downdip_Slip',  pnames))...
%         ))/pi;  % change sign of U2 for negative dip bug in disloc
%     qsig(iq) = NaN;
%     qnames{iq} = sprintf('Derived_Okada2_rake_in_deg_CCW_');
%     
%     % geometric potency
%     iq = iq+1;
%     q0(iq) = p0(get_parameter_index('Okada2_Length',pnames))...
%         *p0(get_parameter_index('Okada2_Width', pnames))...
%         *norm([p0(get_parameter_index('Okada2_RL_Strike_Slip',pnames))...
%         ,         p0(get_parameter_index('Okada2_Downdip_Slip',  pnames))...
%         ,         p0(get_parameter_index('Okada2_Tensile',       pnames))]);
%     q1(iq) = p1(get_parameter_index('Okada2_Length',pnames))...
%         *p1(get_parameter_index('Okada2_Width', pnames))...
%         *norm([p1(get_parameter_index('Okada2_RL_Strike_Slip',pnames))...
%         ,         p1(get_parameter_index('Okada2_Downdip_Slip'  ,pnames))...
%         ,         p1(get_parameter_index('Okada2_Tensile',       pnames))]);
%     qnames{iq} = sprintf('Derived_Okada2_potency_in_m3___');
    %qsig(iq) = NaN;
    %
    %     qsig(iq) = qsig(iq) * psig(get_parameter_index('Okada2_Length',pnames)) ...
    %                 * psig(get_parameter_index('Okada2_Width',pnames)) ;
    % moment, assuming shear modulus = 30 GPa

    % 20130419 Uncertainty is sum of squares
%     qsig(iq) = sqrt((psig(get_parameter_index('Okada2_RL_Strike_Slip'  ,pnames)))^2  ...
%         + (psig(get_parameter_index('Okada2_Downdip_Slip'  , pnames)))^2   ...
%         + (psig(get_parameter_index('Okada2_Tensile', pnames)))^2);
    % 2014013 aboves gives strange results
    qsig(iq) = NaN;
    iq = iq+1;
    if q0(iq-1) > 0; q0(iq) = 3e10*q0(iq-1); else q0(iq) = NaN; end;
    if q1(iq-1) > 0; q1(iq) = 3e10*q1(iq-1); else q1(iq) = NaN; end;
    qnames{iq} = sprintf('Derived_Okada2_moment_in_Nm____');
    qsig(iq) = NaN;
    
    % magnitude
    iq = iq+1;
    if q0(iq-1) > 0; q0(iq) = 2*log10(q0(iq-1))/3 - 6.03; else q0(iq) = NaN; end
    if q1(iq-1) > 0; q1(iq) = 2*log10(q1(iq-1))/3 - 6.03; else q1(iq) = NaN; end
    qnames{iq} = sprintf('Derived_Okada2_Mw______________');
    qsig(iq) = NaN;
end

% geographic coordinates of Mogi source
if abs(psig(get_parameter_index('Mogi1_Volume_Increase_in_m3_____',pnames))) > 0
    [LatMogi10,LonMogi10]=utm2deg(...
        p0(get_parameter_index('Mogi1_Easting_in_m______________',pnames))...
        ,p0(get_parameter_index('Mogi1_Northing_in_m_____________',pnames))...
        ,utmzone0);
    [LatMogi11,LonMogi11]=utm2deg(...
        p1(get_parameter_index('Mogi1_Easting_in_m______________',pnames))...
        ,p1(get_parameter_index('Mogi1_Northing_in_m_____________',pnames))...
        ,utmzone0);
    iq = iq+1;
    q0(iq) = LatMogi10;
    q1(iq) = LatMogi11;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Mogi1_Centroid_latitude__in_deg_');
    iq = iq+1;
    q0(iq) = LonMogi10;
    q1(iq) = LonMogi11;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Mogi1_Centroid_longitude_in_deg_');
end
if abs(psig(get_parameter_index('Mogi2_Volume_Increase_in_m3_____',pnames))) > 0
    [LatMogi20,LonMogi20]=utm2deg(...
        p0(get_parameter_index('Mogi2_Easting_in_m______________',pnames))...
        ,p0(get_parameter_index('Mogi2_Northing_in_m_____________',pnames))...
        ,utmzone0);
    [LatMogi21,LonMogi21]=utm2deg(...
        p1(get_parameter_index('Mogi2_Easting_in_m______________',pnames))...
        ,p1(get_parameter_index('Mogi2_Northing_in_m_____________',pnames))...
        ,utmzone0);
    iq = iq+1;
    q0(iq) = LatMogi20;
    q1(iq) = LatMogi21;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Mogi2_Centroid_latitude__in_deg_');
    iq = iq+1;
    q0(iq) = LonMogi20;
    q1(iq) = LonMogi21;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Mogi2_Centroid_longitude_in_deg_');
end

% Ellipsoid
if isfinite(psig(get_parameter_index('YangPS_semimajor_axis_a__m______',pnames)))==1 ...
        && isfinite(psig(get_parameter_index('YangPS_semiminor_axis_b_in_m____',pnames)))==1
    
    Yang.a  = p0(get_parameter_index('YangPS_semimajor_axis_a__m______',pnames));
    Yang.b  = p0(get_parameter_index('YangPS_semiminor_axis_b_in_m____',pnames));
    Yang.V0 = 4.0 * pi * (Yang.a)^2 * Yang.b;
    
    Yang.a = p1(get_parameter_index('YangPS_semimajor_axis_a__m______',pnames));
    Yang.b = p1(get_parameter_index('YangPS_semiminor_axis_b_in_m____',pnames));
    Yang.x = p1(get_parameter_index('YangPS_Easting_in_m_____________',pnames));
    Yang.y = p1(get_parameter_index('YangPS_Northing_in_m____________',pnames));
    Yang.z = p1(get_parameter_index('YangPS_Depth_in_m_______________',pnames));
    Yang.p = p1(get_parameter_index('YangPS_Excess_Pressure_in_Pa____',pnames));
    Yang.c = p1(get_parameter_index('YangPS_Azimuth_deg_CCW_from_N___',pnames));
    Yang.d = p1(get_parameter_index('YangPS_Plunge_in_degrees________',pnames));
    Yang.V1 = 4.0 * pi * (Yang.a)^2 * Yang.b;
    
    % get UTM coordinates of center
    if isgeo  == 1
        [Yang.latc,Yang.lonc] = utm2deg(Yang.x,Yang.y,utmzone0);
    else
        Yang.latc = NaN;
        Yang.lonc = NaN;
    end
    
    % make prolate spheroid by repeating b axis
    % make depth negative up by negating z
    [Yang.xs,Yang.ys,Yang.zs] = get_ellipsoid(Yang.x,Yang.y,-1.0*Yang.z,Yang.a,Yang.b,Yang.b,Yang.c,Yang.d);
    
    iq=iq+1;
    q0(iq) = Yang.V0;
    q1(iq) = Yang.V1;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('YangPS_Volume_in_cubic_meters___');
    iq=iq+1;
    q0(iq) = NaN;
    q1(iq) = abs(max(max(Yang.zs)));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('YangPS_depth_to_top_in_meters___');
    iq=iq+1;
    q0(iq) = NaN;
    q1(iq) = abs(min(min(Yang.zs)));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('YangPS_depth_to_bottom_in_meters');
    iq=iq+1;
    q0(iq) = NaN;
    q1(iq) = Yang.lonc;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('YangPS_longitude_in_degrees_____');
    iq=iq+1;
    q0(iq) = NaN;
    q1(iq) = Yang.latc;
    qsig(iq) = NaN;
    qnames{iq} = sprintf('YangPS_latitude_in_degrees______');
end

%  PRINT OUT THE DERIVED PARAMETERS
qnames = truncate_parameter_names(qnames);
iq1 = 1;
iq2 = iq;
uqb(iq1:iq2)=NaN;
lqb(iq1:iq2)=NaN;
% 2012 JUL 10
%qsig(iq1:iq2)=NaN;
for i=iq1:iq2
    qflags{i} = 'D#';
    adj = q1(i)-q0(i);
    %    sadj = NaN;
    sadj = abs(adj/qsig(i));
    outfmt = getfmt(q1(i),qnames{i});
    
    fprintf(1        ,outfmt,qflags{i},i+mparam,qnames{i} ,q0(i),q1(i),adj,qsig(i),sadj,(uqb(i)-lqb(i))/2.0);
    fprintf(fidtxtout,outfmt,qflags{i},i+mparam,qnames{i}, q0(i),q1(i),adj,qsig(i),sadj,(uqb(i)-lqb(i))/2.0);
end


PST1.sigma = colvec(psig);
write_pst(PST1,'PST.OUT');

clear h;

save('qsave.mat','iq','iq1','iq2','qflags','qnames','q0','q1','qsig');
%if ismember(ianneal,[1,2])==1
if isfinite(acosts1(1)) == 1
    iok = find(isfinite(acosts1)==1);
    fprintf(1,'Saving %d sets of trial values of parameters and their costs\n',numel(iok));
    acosts1=acosts1(iok);
    trials=trials(iok,:);
    save('trialvals.mat','trials','acosts1','pnames','crit69','cost00','cost0','cost1','runname','p0','p1');
end

save;

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));

return;


