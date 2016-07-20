function ft = time_function(tfunc, ti, tbreaks, metaparams)
%function ft = time_function8(tfunc, ti, tbreaks, metaparams)
% return value of time function f(t) at ti
%
%    inputs:
%          tfunc      == function handle for time function
%          ti         == epoch in years
%          tbreaks    == reference epochs in years
%          metaparams == constants, such as characteristic time for log or exp
%    output
%          ft         == vector containing value of time function
%                        evaluated at each epoch
%
% Metaparameters are pre-defined variables in a temporal function
%    that are not to be solved for in the inversion.  For example, given an
%    exponential model, the reference time epoch and characteristic time
%    constant are metaparams and should be defined by the user.  If your
%    parameterization does not include any metaparams, then metaparams =
%    nan. 
% tbreaks are epochs that define separate intervals in the temporal
%    analysis.  For example, If using a piecewise-lienar parameterization,
%    each epoch should be listed as an element in tbreaks.  If choosing to
%    only model a 2-segment linear model, the first and last epochs as well
%    as 1 epoch in between should be elements in tbreaks.  At minimum (or in
%    other words, if treating the whole time span as a single interval
%    in the parameterization), the first and last epoch should be listed in
%    tbreaks: tbreaks = [min(tu), max(tu)]
%
% Additional parameterizations can be added, including using a combination of parameterizations (see case 'exp2' for an example)
%
% Built-In Parameterizations:
% 'rate'/'secular' - constant rate 
%                      metaparams = nan
%                      tbreaks = [tu(1), tu(end)];
% 'nsegs' - n-segment piecewise linear parameterization (for complete pieceiwse linear, use 'pwl' or 'ber')
%                      metaparams = nan;
%                      tbreaks: dependent on reference epochs
% 'pwl', 'ber' - piecewise linear (breaks at every epoch); piecewise linear with Berardino et al. (2002) parameterization
%                      metaparams = nan; 
%                      tbreaks = tu(1:end);
% 'ber_tikh' - piecewise linear with Berardino et al. (2002) parameterization and Tikhonov regularization
%                      metaparams = [ metaparams(1) = order of regularization;
%                                     metaparams(2) = lower limit for reg. param. (beta) range;
%                                     metaparams(3) = upper limit for reg. param. (beta) range;
%                                     metaparams(4) = delta for reg. param. (beta) range]
%                                   defined such that range of beta = metaparams(1):metaparams(4):metaparams(3); (not used for building temporal function)
%                      tbreaks = tu(1:end)
% 'poly...' - polynomial parameterization of order '...' (for example, 2nd order = 'poly2', 2nd order with constant = 'poly02')
%                      metaparams = nan;
%                      tbreaks = [tu(1), tu(end)]
% 'exprdecay' - exponential decay (rate): f(t) = 1.0-exp(-1*(dt/tau));
%                      metaparams = [ metaparams(1) = tstart; reference epoch
%                                     metaparams(2) = tau; characteristic time constant]
%                      tbreaks = [tu(1), tu(end)];           
% 'exprgrowth' - exponential growth (rate): f(t) 1.0 - exp(1*dt/tau);
%                      metaparams = [ metaparams(1) = tstart; reference epoch
%                                     metaparams(2) = tau; characteristic time constant]
%                      tbreaks = [tu(1), tu(end)]; 
%
% Additional parameterizations can be added, including using a combination of parameterizations (see case 'okmokexp4' for an example relating to test case)
% 
% 2014-JUL-15 Kurt Feigl
% 2015 - 04 - 17 Elena C. Baluyut, UW-Madison
%
% UPDATES:
% 2015-10-12 - add Okmok parameterization for modified exponential
% okmokexp3, ECB 



%    If the persistent variable does not exist the first time you issue
%     the persistent statement, it will be initialized to the empty matrix.

persistent nwarning

% if exist('nwarning','var') == 0
% if nwarning == []
%     nwarning = 0;
% end

if nargin ~= 4
    error(sprintf('wrong number of arguments %d. Need 4\n',nargin));
end

% Start of the time function, before this epoch, time function is 0
tstart=metaparams(1);

% Prune and sort
iok = isfinite(tbreaks);
tbreaks = colvec(sort(unique(tbreaks(iok))));
% Initialize
jj = 0;

%% Sort and find temporal function from tfunc
%   ft represents vector storing time function values
%   jj is looping variable
%   mparams represent number of parameters in parameterization 

if numel(tbreaks) > 0
    switch(lower(tfunc))
        
        case {'rate','secular'} % CONSTANT RATE
            if numel(tbreaks) > 1
                if numel(nwarning) == 0
                    warning('ignoring %d extra tbreaks in secular parameterization',numel(tbreaks));
                    nwarning = 1;
                else
                    nwarning = nwarning + 1;
                end
            end
            mparams = 1; % only one interval
            ft = zeros(mparams,1);
            jj = jj+1;
            ft(jj) = ti - tbreaks(1);
            
        case {'nsegs'} % N SEGMENTS: Piece-wise linear with fewer breaks than epochs
            mparams = numel(tbreaks)-1; % number of intervals
            ft = zeros(mparams,1);
            for j=1:numel(tbreaks)-1
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    %                    ft(jj) = ti - tbreaks(j);
                    if ti < tbreaks(j+1)            % and ends during interval
                        ft(jj) = ti - tbreaks(j);
                    elseif ti >= tbreaks(j+1)        % and ends after interval
                        ft(jj) = tbreaks(j+1) - tbreaks(j);
                    end
                else
                    ft(jj) = 0;
                end
            end
            
        case {'pwl', 'ber', 'ber_tikh'} % PIECEWISE LINEAR: options of Berardino and Tikhonov regularization
            mparams = numel(tbreaks)-1; % number of intervals
            ft = zeros(mparams,1);
            for j=1:numel(tbreaks)-1
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    if ti < tbreaks(j+1)            % and ends during interval
                        ft(jj) = ti - tbreaks(j);
                    elseif ti >= tbreaks(j+1)        % and ends after interval
                        ft(jj) = tbreaks(j+1) - tbreaks(j);
                    end
                else
                    ft(jj) = 0;
                end
            end
            
        case {'poly01','poly02','poly03','poly04','poly05','poly1','poly2','poly3','poly4','poly5', 'poly6', 'poly8'} % POLYNOMIALS
            nintervals = numel(tbreaks)-1; % number of intervals
            norder = str2num(tfunc(5:end));
            if tfunc(5) == '0'  % leading zero includes zero-th order term
                iorder = 0;
                mparams = (norder+1) * nintervals;
            else
                iorder = 1;
                mparams = norder * nintervals;
            end
            fms = time_function_ref('pwl',ti,tbreaks,[]); % call this function to avoid repeating code above
            jj = 0;
            for j=1:numel(fms)
                % zero-th order term
                if iorder == 0
                    jj = jj+1;
                    if abs(fms(j)) > 0.
                        ft(jj) = 1.0;
                    else
                        ft(jj) = 0.0;
                    end
                end
                % higher order terms
                for k = 1:norder
                    jj = jj+1;
                    ft(jj) = power(fms(j),k);
                end
            end
            
        case {'exprdecay'} % exponential growth at exponentially decaying rate
            %time constant set to tstart
            mparams = 1;
            ft = zeros(mparams,1);
            jj = jj+1;
            if ti >= tstart
                dt = ti - tstart;
                ft(jj) = 1.0-exp(-1*dt/metaparams(2)); % exponential decay
            else
                ft(jj) = 0;
            end
            
        case {'exprgrowth'} % exponential with time constant set to tstart
            mparams = 1;
            ft = zeros(mparams,1);
            jj = jj+1;
            if ti >= tstart
                dt = ti - tstart;
                ft(jj) = 1.0 - exp(1*dt/metaparams(2)); % exponential growth
            else
                ft(jj) = 0;
            end
          
        case {'okmokexp3'} % 2 exprdecay with parameterizations with secular rate in between. For Okmok test case. Example of combining time functions.
            jj=2;
            mparams = 2;
            ft = zeros(mparams,1);
            tswitch1=metaparams(3); % switch time
            tswitch2=metaparams(4); % switch time
            if tstart >= tswitch1
                disp('error')
            end
            
            if ti >= tstart && ti < tswitch1
                dt1 = ti - tstart;
                ft(1) = 1.0-exp(-1*dt1/metaparams(2)); % rate decays exponentially
                ft(2) = 0;
            elseif ti >= tswitch1 <= tswitch2
                dt1 = ti - tstart;
                dt2 = ti - tswitch1;
                ft(1) = 1.0-exp(-1*dt1/metaparams(2)); % rate decays exponentially
                ft(2) = dt2; 
            else
                ft(1) = 0;
                ft(2) = 0;
            end    
            
        case {'okmokexp4'} % 2 exprdecay with parameterizations with secular rate in between. For Okmok test case. Example of combining time functions.
            jj=3;
            mparams = 3;
            ft = zeros(mparams,1);
            tswitch1=metaparams(3); % switch time
            tswitch2=metaparams(4); % switch time
            if tstart >= tswitch1
                disp('error')
            end
            
            if ti >= tstart && ti < tswitch1
                dt1 = ti - tstart;
                ft(1) = 1.0-exp(-1*dt1/metaparams(2)); % rate decays exponentially
                ft(2) = 0;
                ft(3) = 0;
            elseif ti >= tswitch1 && ti < tswitch2
                dt1 = ti - tstart;
                dt2 = ti - tswitch1;
                ft(1) = 1.0-exp(-1*dt1/metaparams(2)); % rate decays exponentially
                ft(2) = dt2; % rate
                ft(3) = 0;
            elseif ti >= tswitch2
                dt1 = ti - tstart;
                dt2 = ti - tswitch1;
                dt3 = ti - tswitch2;
                ft(1) = 1.0-exp(-1*dt1/metaparams(2)); % rate decays exponentially
                ft(2) = dt2; % rate
                ft(3) = 1.0-exp(-1*dt3/metaparams(5));%dt3; % rate
            else
                ft(1) = 0;
                ft(2) = 0;
                ft(3) = 0;
            end
            
            
        otherwise
            error(sprintf('undefined tfunc %s',tfunc));
    end
else
    error('Not enough tbreaks');
end

if jj ~= mparams
    error(sprintf('Miscount:  jj (%d) does not equal mparams (%d)',jj,mparams));
end


return
end

