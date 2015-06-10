function [mhat,F,store,energy,count,msig]=anneal5(FUN,bounds,OPTIONS,varargin)
%ANNEAL     [mhat,F,store,energy,count]=anneal5(FUN,bounds,OPTIONS,x1,x2...,xn)
%Simulated annealing algorithm that tries to find a minimum to the function 'FUN'.
%
%INPUTS:
%
%'FUN' specifies the objective function.  This function should accept an input model vector
%as the first input a scalar cost. Additional context-specific argumentscan be passed to the
%objective function by passing them to 'ANNEAL' after 'OPTIONS'.
%
%'bounds' specifies the upper and lower limits that each model parameter can take on.  This matrix
%must have as many rows as model parameters and two columns.
%
%'OPTIONS' specifies a number of annealing options (empty matrix or zeros for defaults):
%
%     OPTIONS(1) = scale of cooling schedule (default = 4).  Higher numbers
%                  produce more exhaustive searches.
%     OPTIONS(2) = number of individual annealing runs (default = 3).  Higher
%                  numbers produce more exhaustive searches and reduce dependency
%                  on correctly guessing critical temperature.
%     OPTIONS(3) = grid spacing (default = 4).  Higher numbers permit finer levels
%                  of parameter discretization.
%     OPTIONS(4) = temperature scale. To tweak this parameter, inspect
%                  a graph of 'energy'.  Values higher than 3 or lower than 1 are
%                  seldom, if ever, warranted.  The default (0) tries several different
%                  values and works well for most problems.
%     OPTIONS(5) = flag that tells the algorithm whether the objective function
%                  can accept matrix input (default = 0); set to '1' if yes.  The
%                  objective function should accept models stored columnwise and
%                  return a vector of costs.  Writing the objective function this
%                  way can increase algorithm speed by 15-20%.
%     OPTIONS(6) = flag that tells the algorithm whether to try improve the solution:
%                  0 = not try to improve (default)
%                  1 = try the Nelder-Mead simplex method using
%                      'fminsearch' function in MATLAB optimization toolbox
%                  2 = try the Constrained optimization method using
%                      'fmincon' function with in MATLAB optimization toolbox
%                       and 'interior-point' algorithm
%     OPTIONS(7) = flag that tells the algorithm whether to display informative
%                  output (default = 0); set to '1' if yes.
%     OPTIONS(8) = flag to use MATLAB's distributed computing routines
%                  set to the number of available processors [default is 1]
%     OPTIONS(9) = flag to Initialize random number generator.
%                  0 == do not initialize [default]
%                  1 == initialize once
%
%OUTPUTS:
%
%'mhat' is the best model found.
%
%'F' is the cost associated with the best model.
%
%'store' is a matrix containing the bestmodels after each sweep.
%
%'energy' is a vector containing the costs corresponding to the models in 'model'.
%
%'count' is the total number of models checked.
%
%Version 1.0  Peter Cervelli 4-26-98.
%Modifications:
% 2007-JUL
%    record results for later use
% 2009-NOV-12 -
%    OPTIONS(8) = 1 to use Matlab distributed computing
% 20100-JUN-21
%    OPTIONS(9) = 1 Initialize random number generator
% 2011-JUN-14
%    adapt to use new optimization
% 2011-DEC-09
%    adapt for Matlab release R2011b
%    try also jacknife

%Check argument syntax

if nargin<2
    help(mfilename)
    %error('Usage: [mhat,F,model,energy,count]=anneal(FUN,bounds,OPTIONS,varargin)')
    error('Incorrect usage')
end

if nargin<3
    OPTIONS=[];
end

if size(bounds,2)~=2
    error('Second argument must be an nx2 matrix of parameter bounds, where n is the number of parameters.');
end


%Check OPTIONS

if isempty(OPTIONS)
    scale=4;
    runs=3;
    grid=4;
    ts=linspace(1.5,2.5,runs);
    matrix=0;
    newton=0;
    talk=1;
    parallel=0;
else
    %OPTIONS(8)=0;
    if OPTIONS(1)
        scale=OPTIONS(1);
    else
        scale=4;
    end
    
    if OPTIONS(2)
        runs=OPTIONS(2);
    else
        runs=3;
    end
    
    if OPTIONS(3)
        grid=OPTIONS(3);
    else
        grid=4;
    end
    
    if OPTIONS(4)
        ts=ones(runs,1)*OPTIONS(4);
    else
        ts=linspace(2,3,runs);
    end
    
    matrix=OPTIONS(5);
    newton=OPTIONS(6);
    if numel(OPTIONS) >= 7
        talk=OPTIONS(7);
    end
    if numel(OPTIONS) >= 8
%         if OPTIONS(8) > 1
%             if talk
%                 fprintf(1,'Checking that matlabpool is open with %d processors....\n',OPTIONS(8));
%             end
%             %matlabpool('open',OPTIONS(8))
%             if matlabpool('size') == OPTIONS(8)
%                 if talk
%                     fprintf(1,'Success.\n');
%                 end
%             else
%                 error(sprintf('Distproc failure.\n'));
%             end
%         end
    end
    if numel(OPTIONS) >= 9
        if OPTIONS(9) == 1
            if talk
                fprintf(1,'%s initializing random number generator.\n',mfilename);
            end
            %        Replace the default stream with a stream whose seed is based on CLOCK, so
            %        RAND will return different values in different MATLAB sessions.  NOTE: It
            %        is usually not desirable to do this more than once per MATLAB session.
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
%             Warning: The RandStream.setDefaultStream static method will be removed in a future release.  Use
%             RandStream.setGlobalStream instead.
        end
    end
end


%Check bounds to make sure they're ok
if max(bounds(:,1)>bounds(:,2))
    error('All the values in the first column of bounds must be less than those in the second.');
end

%Define constants
p=size(bounds,1);   % number of parameters
count=zeros(runs,1);
energy=Inf;
vals=[2.^-(1:grid)];
vals=[vals,0,-vals];
delta=0.5*abs((bounds(:,1)-bounds(:,2)));
kk = 0;
msig = nan(1,p);


% record parameter values if asked
if talk == 2
    %fp = fopen('anneal4.txt','w+t');
    %store = zeros(1,p+2);
    store = nan(p+2,10^5); % arbitrary dimension
end

%Loop through runs
for k=1:runs
    c=0;
    
    if talk
        fprintf(1,'Starting 100 evaluations of objective function to get tc from ts(k) %f\n',ts(k));
    end
    bestmodel=rand(p,100).*((bounds(:,2)-bounds(:,1))*ones(1,100))+bounds(:,1)*ones(1,100);
    
    if matrix
        O=feval(FUN,bestmodel,varargin{:});
    else
        O=zeros(100,1);
        tstart = tic;
        parfor e=1:100
            O(e)=feval(FUN,bestmodel(:,e),varargin{:});
        end
        if talk
            tfor1 = toc(tstart)/100;
            fprintf(1,'Time in seconds per evaluation is          %12.4g\n',tfor1);
            fprintf(1,'Estimated time in seconds to completion is %12.4g\n',2.5e4*tfor1);
        end
    end
    tc=log10(mean(O))-ts(k);
    [v,i]=min(O);
    bestmodel=bestmodel(:,i);
    
    if talk
        fprintf('\n\nBeginning run #%02d. Critical temperature at %12.6g.\n',k,10^tc);
        fprintf('--------------------------------------------------------\n\n');
        fprintf('f-Calls\t\tTemperature\tMinimum f-Value\tDimension NM\n');
        fprintf('--------------------------------------------------------\n');
    end
    
    %Create cooling schedule from critical temperature and scale
    x=scale*[1 2 4 6 10 6 4 2 1];
    t=sum(x);
    temp=logspace(tc+1,tc-1,9);
    T=zeros(t,1);
    C=1;
    for i=1:9
        for j=1:x(i)
            T(C)=temp(i);
            C=C+1;
        end
    end
    
    %Begin annealing
    for w=1:t
        temp=T(w);
        c=c+1;
        
        if talk
            if c/10==floor(c/10)
                fprintf('%7d\t\t%#8.5g\t\t%#8.5g\t%d\n',count(k),temp,min(energy(1:c-1,k)),NM);
            end
        end
        %Visit each parameter
        
        for x=1:p
            if delta(x)
                %Evaluate objective function at each permissible value
                v=bestmodel(x)+vals*delta(x);
                v=v(find((v<=bounds(x,2))&(v>=bounds(x,1))));
                modelmatrix=bestmodel*ones(1,length(v));
                modelmatrix(x,:)=v;
                NM=size(modelmatrix,2);
                count(k)=count(k)+NM;
                %                 if matrix
                %                     O=feval(FUN,modelmatrix,varargin{:});
                %                 else
                O=zeros(NM,1);
                parfor e=1:NM
                    O(e)=feval(FUN,modelmatrix(:,e),varargin{:});
                end
                
                % store values
                %                 if talk == 2
                %                     for e2=1:NM
                %                         fprintf(fp,'%12.6f ',O(e2));
                %                         for k2=1:numel(modelmatrix(:,e2))
                %                             fprintf(fp,'%12.5E ',modelmatrix(k2,e2));
                %                         end
                %                         fprintf(fp,'\n');
                %                     end
                %                 end
                if talk == 2
                    %disp 'dimensions of store:';size(store)
                    for e2=1:NM
                        kk = kk+1;
                        %                         store(kk,1)    = O(e2); % cost
                        %                         store(kk,2)    = temp;
                        %                         %disp 'dimensions of rowvec(modelmatrix(:,e2)):';size(rowvec(modelmatrix(:,e2)))
                        %                         store(kk,3:2+p) = rowvec(modelmatrix(:,e2));
                        store(1,kk)    = O(e2); % cost
                        store(2,kk)    = temp;
                        %disp 'dimensions of rowvec(modelmatrix(:,e2)):';size(rowvec(modelmatrix(:,e2)))
                        store(3:2+p,kk) = modelmatrix(:,e2);
                    end
                end
                
                %Form exponential probability distribution
                [dist,nanflag]=MakePDF(temp,O);
                if nanflag~=0
                    for u=1:length(nanflag)
                        disp(['Warning: the cost function generated a NaN for the following model:']);
                        disp(modelmatrix(nanflag(u)))
                    end
                end
                
                %Sample from probability distribution
                s=find(cumsum(dist)>=rand);
                s=s(1);
                bestmodel(x,1)=modelmatrix(x,s);
                energy(c,k)=O(s);
                model(:,c)=bestmodel;
            end
        end
    end
    
    [F(k,1),i]=min(energy(:,k));
    mhat(:,k)=model(:,i);
    
end

if talk == 2
    %fclose(fp);
    fprintf(1,'Dimensions of stored values %d rows by %d columns\n',size(store,1),size(store,2));
else
    store = model;
end

[F,i]=min(F);
mhat=mhat(:,i);

if exist('msig','var') ~= 1
    msig = nan(size(mhat));
end


% 2011-JUN-25
count = sum(count);

return


function [pdf,nanflag]=MakePDF(temp,v)
%Forms exponential probability distribution given a temperature and vector of costs.
%Internal function for simulated annealing algorithm.

bad=find(isnan(v));
if isempty(bad)
    pdf=eprob(temp,v);
    nanflag=0;
else
    good=find(~isnan(v));
    w=v(good);
    pdf=eprob(temp,w);
    pdf(good)=pdf;
    pdf(bad)=0;
    nanflag=bad;
end
return



function [pdf]=eprob(temp,v)
%Scales cost vector and calculates exponential probability distribution.  The scaling
%is necessary to permit wide ranges in temperature.
%Internal function for simulated annealing algorithm.

toobig=708.3964185322641;
pdf=v/temp;
mpdf=max(pdf);

if mpdf>toobig
    scale=mpdf/toobig;
    pdf=exp(-pdf/scale);
    pdf=pdf/max(pdf);
    pdf=pdf.^scale;
else
    pdf=exp(-pdf);
    pdf=pdf/max(pdf);
end

pdf=pdf/sum(pdf);
return


