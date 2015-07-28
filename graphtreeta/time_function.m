function ft = time_function(tfunc, ti, tbreaks, metaparams)
%function ft = time_function7(tfunc, ti, tbreaks, metaparams)
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
% Additional parameterizations can be added, including using a combination of parameterizations (see case 'exp2' for an example)
%
% 2014-JUL-15 Kurt Feigl
% 2015 - 04 - 17 Elena Baluyut



if nargin ~= 4
    error(sprintf('wrong number of arguments %d. Need 4\n',nargin));
end

% start of the time function, before this epoch, time function is 0
tstart=metaparams(1);

% prune and sort
iok = isfinite(tbreaks);
tbreaks = colvec(sort(unique(tbreaks(iok))));
jj = 0;
if numel(tbreaks) > 0
    switch(lower(tfunc))
        case {'rate','secular'} % CONSTANT RATE
%             if numel(tbreaks) > 1
%                 warning('ignoring extra tbreaks in secular parameterization');
%             end
            mparams = 1; % only one interval
            ft = zeros(mparams,1);
            jj = jj+1;
            ft(jj) = ti - tbreaks(1);
        case {'step','heaviside'} % STEPS
            mparams = numel(tbreaks); % number of breaks
            ft = zeros(mparams,1);
            for j=1:numel(tbreaks)
                jj = jj+1;
                if ti >= tbreaks(j)
                    ft(jj) = 1.0;
                else
                    ft(jj) = 0.0;
                end
            end
        case {'nsegs'} % N SEGMENTS: Piece-wise linear with fewer breaks than epochs
            mparams = numel(tbreaks)-1; % number of intervals
            %mparams = mparams + 1;       % plus 1 for y-intercept ?
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
        case {'bnsegs'} % N SEGMENTS: Piece-wise linear with fewer breaks than epochs
            mparams = numel(tbreaks)-1; % number of intervals
            %mparams = mparams + 1;       % plus 1 for y-intercept ?
            ft = zeros(mparams,1);
            for j=1:numel(tbreaks)-1
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    %                    ft(jj) = ti - tbreaks(j);
                    if ti < tbreaks(j+1) % and ends during interval
                        if tbreaks(j) == metaparams(1) %sets rate between metapaarms(1) and (2) tp zero
                            ft(jj) = 0;
                        else
                            ft(jj) = ti - tbreaks(j);
                        end
                    elseif ti >= tbreaks(j+1)        % and ends after interval
                        ft(jj) = tbreaks(j+1) - tbreaks(j);
                    end
                else
                    ft(jj) = 0;
                end
            end
        case 'step-pwl' % PIECE-WISE LINEAR WITH BREAKS
            nintervals = numel(tbreaks)-1;     % number of intervals
            nbreaks    = numel(tbreaks);       % number of breaks, including first and last
            mparams    = nintervals + nbreaks - 2; % number of parameters
            ft = zeros(mparams,1);
            for j=1:nintervals
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    if ti < tbreaks(j+1)            % and ends during interval
                        ft(jj) = ti - tbreaks(j);
                    elseif ti >= tbreaks(j+1)        % and ends after interval
                        ft(jj) = tbreaks(j+1) - tbreaks(j);
                    end
                else
                    ft(jj) = 0.;
                end
            end
            for j=2:nbreaks-1
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    ft(jj) = 1.;
                else
                    ft(jj) = 0.;
                end
            end
        case 'bradys_zero_expdec'
%             jj=1;
            mparams = 2;
            ft = zeros(mparams,1);
            tswitch=metaparams(2); % last epoch for zero rate
            tstart = metaparams(1);
            %jj = jj+1;
            if ti <= tswitch
            ft(1) = ti - ti;
            ft(2) = 0;
            elseif ti >= tswitch
                dt = ti - tswitch;
                ft(1) = 0;
                ft(2) = 1.0-exp(-1*dt/metaparams(4)); % exponential decay
                
            else
                error('problem with index')
            end
            
        
        case 'step-sec'  % CONSTANT RATE WITH BREAKS
            nintervals = 1;                    % number of intervals
            nbreaks    = numel(tbreaks);       % number of breaks, including first and last
            mparams    = nintervals + nbreaks - 2; % number of parameters
            ft = zeros(mparams,1);
            jj = 1;
            ft(jj) = ti - tbreaks(1);
            for j=2:nbreaks-1
                jj = jj+1;
                if ti >= tbreaks(j)  % interferogram starts during interval
                    ft(jj) = 1.;
                else
                    ft(jj) = 0.;
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
        case {'gaussian'} % Gaussian with standard deviation set to tstart
            %             if numel(tbreaks) > 1
            %                 warning('ignoring extra tbreaks in gaussian parameterization');
            %             end
            mparams = 1; % only one interval
            ft = zeros(mparams,1);
            dt = ti - tbreaks(1);
            jj = jj+1;
            ft(jj) = pdf('normal',dt,0,tstart);
            
          
        case {'sin'}
            tstart = metaparams(1); %reference epoch
            f = 1/metaparams(2); %frequency, in terms of period
            phi = metaparams(3); %shift (if needed)
            A = metaparams(4); %multiplier
            mparams = 1;
            ft = zeros(mparams,1);
            jj = jj+1;
            if ti >= tstart
                dt = ti-tstart;
                ft(jj) = A*sin(2*pi*f*dt+phi);
            else
                ft(jj) = 0;
            end
        case {'cos'}
            tstart = metaparams(1); %reference epoch
            f = 1/metaparams(2); %frequency, in terms of period
            phi = metaparams(3); %shift (if needed)
            mparams = 1;
            ft = zeros(mparams,1);
            jj = jj+1;
            if ti >= tstart
                dt = ti-tstart;
                ft(jj) = cos(2*pi*f*dt+phi);
            else
                ft(jj) = 0;
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
            
        case {'exp2'} % exponential with time constant set to tstart
            jj=2;
            mparams = 2;
            ft = zeros(mparams,1);
            tswitch=metaparams(3); % switch from exponential growth to decay
            if tstart >= tswitch
                disp('error')
            end
            
            if ti < tswitch && ti >= tstart
                dt1 = ti - tstart;
                ft(1) = 1.0 - exp(+1*dt1/metaparams(2)); % exponential growth
                ft(2) = 0;
            elseif ti >= tswitch
                dt1 = tswitch - tstart;
                dt2 = ti - tswitch;
                ft(1) = 1.0 - exp(+1*dt1/metaparams(2)); % exponential growth
                ft(2) = 1.0 - exp(-1*dt2/metaparams(4)); % exponential decay
            else
                ft(1) = 0;
                ft(2) = 0;
            end

        case {'pwl', 'berL'}
            mparams = numel(tbreaks)-1; % number of intervals
            %mparams = mparams + 1;       % plus 1 for y-intercept ?
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
        otherwise
            error(sprintf('undefined tfunc %s',tfunc));
    end
else
    error('Not enough tbreaks');
end

if jj ~= mparams
%    error(sprintf('Miscount:  jj (%d) does not equal mparams (%d)',jj,mparams));
end


return
end

