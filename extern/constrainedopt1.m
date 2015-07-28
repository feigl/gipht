function [mhat,fval,model,energy,count]=constrainedopt1(FUN,bounds,OPTIONS,varargin)
%function [mhat,fval,model,energy,count]=constrainedopt1(FUN,bounds,OPTIONS,varargin)
%Constrained optimization tries to find a minimum to the function 'FUN'.
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
%'model' is a matrix containing the bestmodels after each sweep.
%
%'energy' is a vector containing the costs corresponding to the models in 'model'.
%
%'count' is the total number of models checked.
%
%Version 1.0  Kurt Feigl 2011-JUN-25
% 
% 20150727 FAILS
% Error using sqpLineSearch (line 131)
% Finite difference derivatives at initial point contain Inf or NaN values. Fmincon cannot continue.
% Error in fmincon (line 806)
%     [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = sqpLineSearch(funfcn,X,full(A),full(B),full(Aeq),full(Beq), ...
% Error in constrainedopt1 (line 173)
%     [mstar,fval,exitflag,output]=fmincon(FUN,mhat...
% 
% Error in gipht_step2 (line 529)
%         [p1,f,model,energy,count] = constrainedopt1(objfun,bounds,options,fitfun,DST,PST,TST); 
% 

%Check argument syntax

if nargin<2
    error('Usage: [mhat,F,model,energy,count]=constrainedopt1(FUN,bounds,OPTIONS,varargin)')
end

if nargin<3
    OPTIONS=[];
end

if size(bounds,2)~=2
    error('Second argument must be an nx2 matrix of parameter bounds, where n is the number of parameters.');
end

%Check OPTIONS

if isempty(OPTIONS)
    matrix=0;
    newton=0;
    talk=1;
    parallel=0;
else
    matrix=OPTIONS(5);
    newton=OPTIONS(6);
    if numel(OPTIONS) >= 7
        talk=OPTIONS(7);
    end
    if numel(OPTIONS) >= 8
        if OPTIONS(8) > 1
            fprintf(1,'Checking that matlabpool is open with %d processors....\n',OPTIONS(8));
            %matlabpool('open',OPTIONS(8))
            if matlabpool('size') == OPTIONS(8)
                fprintf(1,'Success.\n');
            else
                error(sprintf('Distproc failure.\n'));
            end
        end
    end
end

%Check bounds to make sure they're ok

if max(bounds(:,1)>bounds(:,2))
    error('All the values in the first column of bounds must be less than those in the second.');
end

%Define constants
p=size(bounds,1);
energy=Inf;


mhat = (bounds(:,1)+bounds(:,2))/2.0;
model = zeros(size(bounds(:,1)));

% find fixed parameters
mparams = numel(bounds)/2
delta=abs((bounds(:,1)-bounds(:,2)));
ifixed = find(delta < eps);
Aeq = zeros(mparams,mparams);
Aeq(ifixed,ifixed) = 1;
beq = zeros(size(mhat));
beq(ifixed) = mhat(ifixed);

count = 0;
cost0=feval(FUN,mhat,varargin{:})
lb = bounds(:,1);
ub = bounds(:,2);



% if newton
%     if newton==1
%         if numel(which('fminsearch')) > 0
%             if talk
%                 fprintf('\nStarting Simplex method FMINSEARCH.\n\n')
%             end
%             % Use old routine
%             %mstar=fmins(FUN,mhat(:,k),[],[],varargin{:});
%             % Use new routine
%             %objfun,bounds,options,fitfun,DST,PST,TST
%             %function [x,fval,exitflag,output] = fminsearch(funfcn,x,options,varargin)
%             %             Exiting: Maximum number of function evaluations has been exceeded
%             %          - increase MaxFunEvals option.
%             %          Current function value: 0.106743
%             
%             [mstar,fval,exitflag,output]=fminsearch(FUN,mhat(:,k),optimset('MaxFunEvals',1e5),varargin{:});
%             if exitflag == 1 || (exitflag == 0 && fval < F(k,1))
%                 cost1=feval(FUN,mstar,varargin{:});
%                 if cost1<F(k,1)
%                     if mstar>bounds(:,1) & mstar<bounds(:,2)
%                         mhat(:,k)=mstar;
%                         F(k,1)=cost1;
%                         if talk
%                             fprintf('\nSimplex method FMINSEARCH lowered cost to %10.4g and remained within constraints.\n\n',cost1)
%                         end
%                     else
%                         if talk
%                             fprintf('\nSimplex method FMINSEARCH lowered cost but failed to remain within constraints.\n\n')
%                         end
%                     end
%                 else
%                     if talk
%                         output
%                         fprintf('\nSimplex method FMINSEARCH failed with exitflag =  %d.\n\n',exitflag)
%                     end
%                 end
%             end
%         else
%             warning('missing FMINCON routine from MATLAB optimization toolbox')
%         end
%     elseif newton==2
if numel(which('fmincon')) > 0
    if talk
        fprintf('\nStarting Constrained method FMINCON.\n\n')
    end
    %Old Routine
    %mstar=constr(FUN,mhat(:,k),[],bounds(:,1),bounds(:,2),[],varargin{:});
    % New Routine
    %             x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    %             x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    %             x = fmincon(problem)
    %             [x,fval] = fmincon(...)
    %             [x,fval,exitflag] = fmincon(...)
    %             [x,fval,exitflag,output] = fmincon(...)
    
   % optoptions = optimset('Algorithm','interior-point','Display','Iter'); % run interior-point algorithm
  optoptions = optimset('Algorithm','sqp','Display','Iter'); % run interior-point algorithm

    [mstar,fval,exitflag,output]=fmincon(FUN,mhat...
        ,[],[],Aeq,beq,lb,ub...
        ,[],optoptions,varargin{:});
    if exitflag == 1 || exitflag == 0
        cost1=feval(FUN,mstar,varargin{:});
        if cost1<cost0
            mhat = mstar;
            fval = cost1;
            if talk
                fprintf('\nConstrained optimization lowered cost to to %10.4g \n\n',cost1)
            end
        end
    else
        fval = cost0;
        if talk
            output
            fprintf('\nConstrained method FMINCON failed with exitflag =  %d.\n\n',exitflag)
        end
    end
else
    warning('missing FMINCON routine from MATLAB optimization toolbox')
end

return


