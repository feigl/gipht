function [TSTP, PST1, trials] = anneal_with_surrogate(objfun, DST, PST, TST, maxiter, verbose, imode, nf, runname, printfun, options)

fobjval = nan(maxiter,1);
fobjval(1) = 1.e99;
for iter = 1:maxiter
    fprintf(1,'Beginning iteration number %3d of %3d\n',iter,maxiter);
    %         % read or generate matrix of partial derivatives
    fnametst = 'TSTP.DAT';
    %         if fexist(fnametst)
    %             TSTP = read_tst(fnametst);
    %             [ndum,mdum] = size(TSTP.partial_wrt_1param);
    %             if [ndum,mdum] == [ndata, mparam]
    %                 igenptl = 0;
    %             else
    %                 warning(sprintf('Existing file named %s has wrong size. %d rows %d columns \n',ndum,mdum));
    %                 igenptl = 1;
    %             end
    %         else
    %             igenptl = 1;
    %         end
    
    % generate partial derivatives
    fprintf(1,'Generating partial derivatives for 1st-order Taylor expansion\n');
    % % make scale match difference between bounds
    %   PST.scale = colvec(ub-lb)/100.;
    PST.scale = (PST.ub - PST.lb)/10.;
    
%     TST
    PST.fitfun
    [rngdum,TSTP] = generate_partials(DST,PST,TST);
%     TSTP
%     whos TSTP.partial_wrt_1param
%     size(TSTP.partial_wrt_1param)
    
    %             % write to a file and check
    %             ierr = write_tst(TSTP,fitfun,fnametst);
    %             TSTP2 = read_tst(fnametst);
    %             nerr = numel(find(abs(TSTP.partial_wrt_1param-TSTP2.partial_wrt_1param) > 1.e-6));
    %             if nerr > 0
    %                 warning(sprintf('Found %d large numerical errors in partial derivatives\n',nerr));
    %             end
    
    % check for NaNs
    [ibad,jbad] = find(isfinite(TSTP.partial_wrt_1param) == 0);
    if numel(ibad) > 0
        warning(sprintf('Found %d NaN values in matrix of partial derivatives. Replacing with zeros.\n',numel(ibad)));
        %     for j=1:numel(ibad)
        %         fprintf(1,'ibad = %10d jbad = %10d\n',ibad(j),jbad(j));
        %     end
        TSTP.partial_wrt_1param(ibad,jbad) = 0;
    end
    
    % evaluate at randomly perturbed estimate of parameters  within bounds
    if verbose == 2
        PST1 = PST;
        %PST1.p1 = PST1.p0 + (rand(size(PST1.p0))-0.5*ones(size(PST1.p0))) .* (ub-lb);
        PST1.p1 = PST1.p0 + (rand(size(PST1.p0))-0.5*ones(size(PST1.p0))) .* (PST.ub-PST.lb);
        %         for i=1:mparam
        %            fprintf(1,'%3d %s %10.3g\n',i,char(pnames{i}),PST1.p1(i)-PST1.p0(i));
        %         end
        % evaluate exactly
        mdl1e = feval(PST.fitfun,DST,PST1,TST);
    end
    %
    % now use the Taylor Expansion to the Fitting function
    %fitfun_orig = fitfun;
    %fitfun = 'funtaylor1';
    %     PST.fitfun = 'funtaylor1';
%    TST = TSTP;
%     clear TSTP;
    
    if verbose == 2
        % now evaluate approximately
        %mdl1p = feval('funtaylor1',DST,PST1,TST);
        mdl1p = feval('funtaylor1',DST,PST1,TSTP);
        %clear PST1;
        
        epdiff = nanstd(mdl1p-mdl1e)/2./pi;
        fprintf(1,'Standard deviation of difference between exact and approximate: %.5f [cycles]\n',epdiff);
        % figure comparing Exact vs. Aprox
        nf=nf+1;h(nf)=figure;subplot(2,1,1);
        plot(mdl1e/2./pi,'b-'); hold on;
        plot(mdl1p/2./pi,'r-');
        ylabel('model phase value (cycles)');
        legend('exact fitting function','1st order Taylor expansion','location','best');
        subplot(2,1,2);
        plot((mdl1p-mdl1e)/2./pi,'g-');
        xlabel('index');ylabel('phase value (cycles)');legend('aprox error','location','best');
        title(sprintf('standard deviation = %.4f cycles\n',epdiff));
        feval(printfun,sprintf('%s_Taylor1Check',runname));
        clear mdl1p, mdl1e;
    end
    
    % perform annealing if required
    if imode >= 2
        
        % set bounds
        bounds(:,1) = PST.lb;
        bounds(:,2) = PST.ub;

        fprintf (1,'\nStarting anneal5 without recording...\n');
        options(7) = 1;acosts1(1)=NaN;
        tanneal=tic;
        %[p1,f,model,energy,count] = anneal5(objfun,bounds,options,fitfun,DST,PST,TST);
        %[p1,f,trials,energy,count] = anneal5(objfun,bounds,options,fitfun,DST,PST,TST);
        [p1,f,trials,energy,count] = anneal5(objfun,bounds,options,@funtaylor1,DST,PST,TSTP);

        fprintf (1,'\nAnnealing ended after %15.0f seconds\n',toc(tanneal));
        
        % record value of objective function
        fobjval(iter) = f
        
        % should do bootstrap here
        msig = nan(size(p1));
        
        % update estimate if solution is better
        if f < nanmin(fobjval)
            fprintf(1,'FOUND a better solution with objective function = %g\n',f);
            PST.p0   = colvec(p1);
        end
    end
    
    
end % loop over iterations

PST1 = PST;

if imode >= 2
    % save final estimate
    PST1.p1    = colvec(p1);
    PST1.sigma = colvec(msig);
else
    trials = nan;
end

return
