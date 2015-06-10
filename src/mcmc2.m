function [mout,mMAP,accrate]=mcmc2(fitfun,niter,mstart,msig,lb,ub,xdat,yobs,ysig)
%function [mout,mMAP,accrate]=mcmc2(fitfun,niter,mstart,msig,lb,ub,xdat,yobs,ysig)
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
%
% The main loop.
fprintf(1,'Starting %s\n',mfilename);
tstart = tic;
fprintf(1,'%12s %12s %12s %12s\n','Iteration','Naccepted','logMAP','accrate');
for k=2:niter
    % Generate a candidate model from the previous model.
    mcandidate=feval('generate2',mcurrent,msig);
  
    % Evalate the logarithm of the acceptance ratio.
    lpcandidate=feval('logprior2',mcandidate,lb,ub);
    llcandidate=feval('loglikelihood2',mcandidate,xdat,yobs,ysig,fitfun);
    %Next two lines give equal results because standard normal distribution is symmetric about zero
    lr1=feval('logproposal2',mcandidate,mcurrent,msig);
    lr2=feval('logproposal2',mcurrent,mcandidate,msig);
    lpcurrent=feval('logprior2',mcurrent,lb,ub);
    llcurrent=feval('loglikelihood2',mcurrent,xdat,yobs,ysig,fitfun);
    % not symmertic 
    logalpha=lpcandidate+llcandidate+lr1-lpcurrent-llcurrent-lr2;
    % can simplify for normal distribution because it is symmetric
    %logalpha=lpcandidate+llcandidate     -lpcurrent-llcurrent;
    %
    % Take the minimum of the log(alpha) and 0.
    %
    if (logalpha>0)
        logalpha=0;
    end
    
    %
    % Generate a U(0,1) random number and take its logarithm.
    %
    logt=log(rand());
    %
    % Accept or reject the step.
    %
    if (logt < logalpha)
        %
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
        %
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
fprintf(1,'Leaving %s after %10.4e seconds\n',mfilename,toc(tstart));

return
end
