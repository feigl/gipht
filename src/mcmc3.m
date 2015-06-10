function [PSTout,h,mskip] = mcmc3(PSTin,DST,TST,mstepscale,nburnin,nskip,niter)
% PST == Parameter structure
mstart = PSTin.p0;   % initial estimate
lb     = PSTin.lb;   % lower bound
ub     = PSTin.ub;   % upper bound
mstep  = mstepscale * (ub-lb); % step size for proposals
fitfun = PSTin.fitfun; % fitting function
mstart = PSTin.p0;
% DST == Data structure

%function [mout,mMAP,accrate]=mcmc2(fitfun,niter,mstart,mstep,lb,ub,xdat,yobs,ysig)
% modified 20140107 by Kurt Feigl
% from
% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
%
% mout=mcmc(logprior,loglikelihood,generate,logproposal,m0,niter)
%
%  logprior           Name of a function that computes the log of
%                     the prior distribution.
%  loglikelihood      Name of a function the computes the log of
%                     the likelihood.
%  generate           Name of a function that generates a random
%                     model from the mcurrent model using the
%                     proposal distribution.
%  logproposal        Name of a function that computes the log of
%                     the proposal distribution r(x,y).
%  m0                 Initial model.
%  niter              Number of iterations to perform.
%
%  mout               MCMC samples.
%  mMAP               Best model found in the MCMC simulation.
%  accrate            Acceptance rate
%
%
%  Figure out some size information.
%
mparams=length(mstart);
%
% Allocate space for the results.
%
mout=zeros(mparams,niter);
%
% Initialize the starting point.
%
mout(:,1)=mstart;
mcurrent=mstart;
lMAP=-Inf;
mMAP=mcurrent;
nacc=0;
accrate=0;

matfilename = sprintf('%s_%s.mat',mfilename,fitfun);

for i=1:mparams
    db = (ub(i)-lb(i))/2.;
    if db > 0.
        fprintf(1,'%03d %32s %12.4g +/- %12.4g\n',i,char(PSTin.names{i}),mstart(i),db);
    end
end

%
% The main loop.

if fexist(matfilename) ~= 1
    fprintf(1,'Starting %s for %d iterations\n',mfilename,niter);
    tstart = tic;
    for k=2:niter
        
        if mod(k,1000) == 0
            %fprintf(1,'%12d %12d %#12.4g %#12.4g %#12.4g\n',k,nacc,lMAP,logalpha,accrate);
        end
        if k == 100
            secperit = toc(tstart)/k;
            tremain = (niter-k)*secperit;
            fprintf(1,'Time in seconds per iteration = %12.4f\n',secperit);
            fprintf(1,'Time in seconds remaining     = %12.4f\n',tremain);  
            fprintf(1,'%12s %12s %12s %12s %12s\n','Iteration','Naccepted','lMAP','accrate','sec per it');
        end
        % Generate a candidate model from the previous model.
        mcandidate=feval('generate2',mcurrent,mstep);
        
        % Evalate the logarithm of the acceptance ratio.
        lpcandidate=feval('logprior2',mcandidate,lb,ub);
        %llcandidate=feval('loglikelihood2',mcandidate,xdat,yobs,ysig,fitfun);
        PSTin.p0 = mcandidate;
        llcandidate=feval('loglikelihood3',mcandidate,fitfun,DST,PSTin,TST);
        %Next two lines give equal results because standard normal distribution is symmetric about zero
        %lr1=feval('logproposal2',mcandidate,mcurrent,mstep)
        %lr2=feval('logproposal2',mcurrent,mcandidate,mstep)
        lpcurrent=feval('logprior2',mcurrent,lb,ub);
        %llcurrent=feval('loglikelihood2',mcurrent,xdat,yobs,ysig,fitfun);
        llcurrent=feval('loglikelihood3',mcurrent,fitfun,DST,PSTin,TST);
        % not symmertic
        %logalpha=lpcandidate+llcandidate+lr1-lpcurrent-llcurrent-lr2
        % can simplify for normal distribution because it is symmetric
        logalpha=lpcandidate+llcandidate     -lpcurrent-llcurrent;
        %
        % Take the minimum of the log(alpha) and 0.
        %
%         if (logalpha>0)
%             logalpha=0;
%         end
        if isnan(logalpha) == 1
            logalpha = 0;
        end
        
        %
        % Generate a U(0,1) random number and take its logarithm.
        %
        logt=log(rand());
        %
        % Accept or reject the step.
        if (logt < logalpha)
            %fprintf(1,'%12d %12d %#12.4g %#12.4g %#12.4g\n',k,nacc,lMAP,accrate,toc(tstart)/k);         

            % Accept the step.
            %
            mcurrent=mcandidate;
            nacc=nacc+1;
            %
            % Update the MAP solution if this one is better.
            %
            if ((lpcandidate+llcandidate) > lMAP)
                lMAP=lpcandidate+llcandidate;
                mMAP=mcandidate;
                % Tell us something.
                fprintf(1,'%12d %12d %#12.4g %#12.4g\n',k,nacc,lMAP,accrate);
               end
        else
            %
            % Reject the step.
            if mod(k,1000) == 0
                fprintf(1,'%12d %12d %#12.4g %#12.4g %#12.4g\n',k,nacc,lMAP,accrate,toc(tstart)/k);
            end         
        end
        %
        % Record the result.
        %
        mout(:,k)=mcurrent;
        accrate=nacc/k;
        %     if accrate >= 0.5
        %         break
        %     end
    end
    fprintf(1,'%12d %12d %12.4f %12.4f\n',k,nacc,lMAP,accrate);
    
    fprintf(1,'Completed %d iterations after %10.4e seconds\n',k,toc(tstart));
    
    
    save(matfilename);
end

load(matfilename);

%downsample results to reduce correlation
k=(nburnin:nskip:niter);

mskip=mout(:,k);
[ndumb,nsamp] = size(mskip);

%histogram results, and find the modes of the subdistributions as an
%estimte of the MAP model
% disp(['m_map','  m_true'])
% [mMAP,mtrue]

% estimate the credible intervals
alpha = 2.0 * cdf('normal',-1,0,1); % 68 percent confidence
%alpha = 2.0 * cdf('normal',-2,0,1); % 95 percent confidence
cred = 1.0 - alpha;
mbot = nan(mparams,1);
mtop = nan(mparams,1);
msig = nan(mparams,1);
mmid = nan(mparams,1);
fprintf(1,'%.0f percent credible interval: center  +/- \n',100*cred);

for i=1:mparams
    msort=sort(mskip(i,:));
    mbot(i)  =  msort(round(      (alpha/2.0)*length(mskip)));
    mtop(i)  =  msort(round((1.0 - alpha/2.0)*length(mskip)));
    mmid(i)  =  msort(round(              0.5*length(mskip)));
    msig(i)  =  (mtop(i)-mbot(i))/2.0;
    %fprintf(1,'%.0f percent credible interval for parameter %d is %12.4g to %12.4g\n',100*cred,i,mbot(i),mtop(i));
    fprintf(1,'%03d %32s %12.4g +/- %12.4g\n',i,char(PSTin.names{i}),mmid(i),msig(i));
end



% write estimates for return
PSTout = PSTin;
PSTout.p1 = colvec(mmid);
PSTout.sigma = colvec(msig);


% find indices of free parameters 
jfree= find((ub-lb) > 0.);
mskipfree=mskip(jfree,:);
for i=1:numel(jfree)
   freenames{i} = PSTin.names{jfree(i)};
end

% Draw figures
nf=0;

% %plot a scatter plot and histogram of the posterior distribution

nf=nf+1;h(nf)=figure;
plotmatrix(mskipfree');
printpdf(sprintf('%s_%s_scatter1.pdf',mfilename,fitfun));

% figure;
% gplotmatrix(mskipfree',[],[1:nsamp],'b',[],[],'off','hist',freenames);
% printpdf(sprintf('%s_%s_scatter2.pdf',mfilename,fitfun));

%
%
% %plot parameter correlations
laglen=50;
lags=(-laglen:laglen)';
%for i=1:mparams
for j=1:numel(jfree)
    i = jfree(j);
    nf=nf+1;h(nf)=figure;
    acorr(:,i)=calc_corr(mskip(i,:)',laglen);
    subplot(2,1,1);
    bookfonts;
    plot([0 laglen],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',3);
    hold on
    plot(lags(laglen+1:laglen*2+1),acorr(laglen+1:laglen*2+1,i),'ko');
    hold off
    ylabel(['A ( m_',num2str(i),')']);
    xlabel('Lag');
    ylabel('Corr');
    title(sprintf('parameter %d %s',i,char(PSTin.names{i})),'Interpreter','None');
    
    subplot(2,1,2);
    bookfonts;
    % histogram of all results
    %hist(mout(i,:));
    % histogram of downsampled results, skipping
    hist(mskip(i,:));
    axis([lb(i) ub(i) -10 +Inf]);
    hold on;
    plot([mbot(i) mtop(i)],[-10 -10],'r-','LineWidth',10);
    %plot(mtrue(i),-10,'o','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k');
    plot(mMAP(i) ,-10,'s','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k');
    xlabel(sprintf('%s',char(PSTin.names{i})),'Interpreter','None');
    ylabel('n');
    title(sprintf('%12.4g +/- %12.4g [%12.4g to %12.4g]\n',mmid(i),msig(i),mbot(i),mtop(i)));
    printpdf(sprintf('%s_%s_param%03d.pdf',mfilename,fitfun,i));
end

% %examine the likelihood probability distribution
% figure
% bookfonts
% for i=1:length(mskip)
%     logli(i)=loglikelihood2(mskip(:,i),xdat,yobs,ysig,fitfun)+logprior2(mskip(:,i),lb,ub);
% end
% hist(logli,50);
% % h = findobj(gca,'Type','patch');
% % set(h,'FaceColor','k')
% hold on
% mMAP_li=loglikelihood2(mMAP,xdat,yobs,ysig,fitfun)+logprior2(mMAP,lb,ub);
% %mtrue_li=loglikelihood2(mtrue,xdat,yobs,ysig,fitfun)+logprior2(mtrue,lb,ub);
% plot(mtrue_li,0.,'o','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k');
% %plot(mMAP_li, 0.,'s','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k');
% %legend('MCMC','true','MAP','location','NorthWest');
% xlabel('ln(Likelihood)')
% ylabel('n')
% %xlim([-15 -7.9])
% title('histogram of the log of the likelihood of each sample');
% printpdf(sprintf('%s_%s_histloglike.pdf',mfilename,fitfun));



return
end
