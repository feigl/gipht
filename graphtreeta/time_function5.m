function ft = time_function5(tfunc, ti, tbreaks, metaparams)
%function ft = time_function5(tfunc, ti, tbreaks, metaparams)
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
% 2014-JUL-07 Kurt Feigl
%



if nargin ~= 4
    error(sprintf('wrong number of arguments %d. Need 4\n',nargin));
end

% prune and sort
iok = isfinite(tbreaks);
tbreaks = colvec(sort(unique(tbreaks(iok))));
jj = 0;
if numel(tbreaks) > 0
    switch(lower(tfunc))
        case {'rate','secular'} % CONSTANT RATE
            if numel(tbreaks) > 1
                warning('ignoring extra tbreaks in secular parameterization');
            end
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
        case {'pwl'} % PIECEWISE LINEAR
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
         case {'poly01','poly02','poly03','poly04','poly05','poly1','poly2','poly3','poly4','poly5'} % POLYNOMIALS 
            nintervals = numel(tbreaks)-1; % number of intervals
            norder = str2num(tfunc(5:end));
            if tfunc(5) == '0'  % leading zero includes zero-th order term
                iorder = 0;
                mparams = (norder+1) * nintervals;
            else
                iorder = 1;
                mparams = norder * nintervals;
            end
            fms = time_function5('pwl',ti,tbreaks,[]); % call this function to avoid repeating code above
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
         case {'gaussian'} % Gaussian with standard deviation set to metaparams(1)
%             if numel(tbreaks) > 1
%                 warning('ignoring extra tbreaks in gaussian parameterization');
%             end
            mparams = 1; % only one interval
            ft = zeros(mparams,1);
            dt = ti - tbreaks(1);
            jj = jj+1;
            ft(jj) = pdf('normal',dt,0,metaparams(1));
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

