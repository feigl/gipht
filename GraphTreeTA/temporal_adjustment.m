function [pest, psig, mse, dmod, tfit, pfit, sfitl, sfitu, rd, V, G, sswr, Vx, var, res_n, Ginv] = temporal_adjustment(data,data_sigma,tm,ts,tbreaks,tfunc,metaparams,Ginv,verbose)
%function [pest, psig, mse, dmod, tfit, pfit, sfitl, sfitu, rd, V, G, sswr, Vx, var, res_n] = temporal_adjustment(data,data_sigma,tm,ts,tbreaks,tfunc,metaparams,verbose)
% Given time tags and data values, perform temporal adjustment. Solve the weighted least squares
% problem Gm = d, where G is the design matrix, d is the data vector
% (Qdiff), and m is the unknown parameter vector.  When m is not
% displacement per epoch, inversion results are interpolated to find
% displacement per epoch. Graph theory is used to define the full data
% covariance matrix, whose inverse is used as a weighting matrix.
%
% Avaliable Models:
%     Piecewise-linear ('pwl'): If 1 distinct component, solves directly for
%       diplacement per epoch (setting a zero displacement constraint at the
%       first epoch).  Otherwise, defaults to Berardino et al. (2002)
%       parameterization.
%     Piecewise-linear with Berardino et al. (2002) method ('ber'): solves for
%       vector of rates and interpolates to get displacement per epoch
%     Piecewise-linear with Berardino et al. (2002) method and Tikhonov
%      regularization ('ber_tikh'): same parameterization as 'ber' but adds
%       zero-, first-, or second- order Tikhonov regularization.  Choice to
%       include a weighting to the design matrix
%     Overdetermined temporal parameterizations: included in a library of
%       temporal functions (time_function.m) but can also be defined by user.
%       Must be overdetermined.
%       Solve inversion for coefficients of temporal function and interpolate
%       to find displacement per epoch.
%
%
% INPUTS:
%     data = pairwise differential data for analysis (e.g., Qdiff)
%     data_sigma = standard deviations of pairwise data (e.g., Qdsig)
%     tm = master epochs (chronologically first epochs in pairs)
%     ts = slave epochs (chronologically second epochs in pairs)
%     tbreaks = user defined breaks (in years) for temporal model
%     tfunc = temporal function
%     metaparams = pre-defined parameters corresponding to temporal function
%
% OUTPUTS:
%     pest = estimated least squares solution of unknown parameters
%     psig = estimated uncertainties of solution parameters
%     mse = mean squared error of residuals
%     dmod = modeled value of differential quantity (from 1st epoch (tm) to 2nd epoch (ts))
%     tfit = vector of time values for modeled fit (temporal function)
%     pfit = vector of modeled (epoch-wise) differential values for modeled fit (temporal function)
%     sfitl = lower 1-sigma bounds for modeled fit (temporal function)
%     sfitu = upper 1-sigma bounds for modeled fit (temporal function)
%     rd = rank deficiency
%     V = covariance matrix of pairwise data
%     G = design matrix
%     sswr = sum of squared weighted residuals
%     Vx = covariance of estimated parameters
%     var = sample variance
%
% Elena Baluyut and Kurt Feigl University of Wisconsin-Madison
% 2020/04/24 add fast option to return Ginv

%% Initialize


narginchk(7,9);
nargoutchk(14,16);

% Set metaparams if not defined
if exist('metaparams','var') == 0
    metaparams = nan;
end

% Set verbose if not defined
if exist('verbose','var') == 0
    verbose = 1;
end

% Set dofast if not defined
if exist('Ginv','var') == 1
    dofast = 1;
else
    dofast = 0;
end

if verbose == 1
    fprintf(1,'%s begins ...\n',mfilename);
end

% Check for proper inversion function
iok=which('ls_with_cov');
if numel(iok) < 1
    error('Need ls_with_cov function')
end

% Find number of data and define variables
%ndat = numel(data);
[ndat,ncols] = size(data);
if ncols ~= 1
    error('Data vector has more than 1 column\n');
end
if ndat < 1
    error('Number of data must be at least 1.\n');
end
if numel(data_sigma) == ndat
    dsig = colvec(data_sigma);
else
    error('data_sigma has wrong dimensions');
end
dmod = nan(size(data)); % set vector used for storing modeled displacement

if dofast == 1
    if verbose == 1
        fprintf(1,'Starting fast option.\n');
    end
end

% make column vectors of epochs
tm = colvec(tm); % master epoch
ts = colvec(ts); % slave epoch
tu = colvec(sort(unique([tm; ts]))); % unique set of epochs
me = numel(tu);  % number of epochs

if verbose
    fprintf(1,'Number of unique epochs me = %d\n',me);
end

if  nanmin(tbreaks) <  nanmin(tu)
    error('min(tbreaks) before first observation');
end

if max(tbreaks) > max(tu)
    max(tbreaks)
    max(tu)
    error('max(tbreaks) after last observation');
end

%% find components
[trees, Q] = findtrees(tm,ts,verbose); % find components and edge-vertex incidence matrix
[ntrees,ndummy] = size(trees); % find number of components
pairs = pairlist(Q); % store indices of epochs in pairs

% % Save the original incidence matrix, to be used later to create
% % correlation and covariance matrices
% Q0 = Q;
% dsig0 = dsig;
% data0 = data;

% Sort reference epochs
tbreaks = colvec(sort(unique(tbreaks)));

%Print reference epochs
if verbose == 1
    fprintf(1,'i,tbreaks(i)\n');
    for i=1:numel(tbreaks)
        fprintf(1,'%d %12.5f\n',i,tbreaks(i));
    end
end

%% Check number of components for PWL
if strcmp(tfunc, 'pwl') == 1 && ntrees > 1
    error('More than one distinct component (rank deficiency greater than 2.\n Consider changing tfunc to rate parameterization outlined by Berardino et al. (2002)\n')
end

%% Build design matrix G

if strcmp(tfunc,'pwl') == 1 && ntrees == 1 && dofast == 0
    % In this case, we can use the incidence matrix Q as the desgin
    % matrix G with an additional constraint of first epoch at 0
    % displacement
    [ndummy,mparams] = size(Q);
    G = Q;
    %set constraint that first epoch has zero displacement
    G(end+1, 1) = 1;
    data(end+1) = 0;
    dsig(end+1) = 1; % aribtrarily assign uncertainty for constraint
elseif strcmp(tfunc, 'ber') == 1  && dofast == 0
    [ G ] = ber_mat( Q, tu);
    mparams = me-1;
else
    switch tfunc
        case {'nseg0','forge'} % need t2nd for these time functions
            t2nd = nanmin([tm; ts]);
        otherwise
            t2nd = nan;
    end
    %For all other cases, use temporal adjustment
    % get number of parameters from temporal function
    %mparams = numel(time_function(tfunc, min([tm; tu]), tbreaks, metaparams));
    mparams = numel(time_function(tfunc, nanmin(tu), tbreaks, metaparams, t2nd));
    if verbose == 1
        fprintf(1,'number of parameters mparams = %d\n',mparams);
    end
    
    % Set up design matrix G
    G = zeros(ndat,mparams);
    
    % Check that system is overdetermined (unless Berardino case)
    if (strncmp(tfunc, 'pwl', 3) == 0) && (strncmp(tfunc, 'ber', 3) == 0)
        if ndat <= mparams
            error('Underdetermined system; cannot perform temporal adjustment.  Change parameterizations.')
        end
    end
    
    % Build design matrix G from temporal function
    switch tfunc
        case {'nseg0','forge'} % need t2nd for these time functions
            for i=1:ndat
                t2nd = ts(i);
                tfm = time_function(tfunc, tm(i), tbreaks, metaparams,t2nd); % time function for master
                tfs = time_function(tfunc, ts(i), tbreaks, metaparams,t2nd); % time function for slave
                for j=1:mparams
                    G(i,j) = tfs(j)-tfm(j); % time function for pair
                end
            end % loop over data
        otherwise
            for i=1:ndat
                tfm = time_function(tfunc, tm(i), tbreaks, metaparams); % time function for master
                tfs = time_function(tfunc, ts(i), tbreaks, metaparams); % time function for slave
                for j=1:mparams
                    G(i,j) = tfs(j)-tfm(j); % time function for pair
                end
            end % loop over data
    end
end

%% Print Data Summary
if verbose == 1
    fprintf(1,'i,  data(i),  dsig(i)\n');
    for i = 1:ndat
        fprintf (1,'%5d %12.4f %12.4f \n',i,data(i),dsig(i));
    end
end

nm = size(G);
rd = mparams - rank(G);
rankG  = rank(G);
if verbose == 1
    fprintf(1,'Length of data vector, excluding constraints %d\n',ndat);
    fprintf(1,'Dimensons of design matrix, including constraints: %d by %d\n',nm(1),nm(2));
    fprintf(1,'Rank of G = %d\n',rankG);
    fprintf(1,'Rank deficiency, including constraints %d\n',rd);
end

% if verbose == 1
%     G
%     figure;
%     spy(G);
%     xlabel('column');ylabel('row');
%     title(sprintf('G matrix for %s has ndat = %d rows and mparams = %d columns'...
%         ,strrep(mfilename,'_',' '),ndat,mparams));
% end


if rd > 0
    if ntrees > 1
        warning(sprintf('Number of trees (%d) is greater than 1.\n',ntrees));
        if strcmp(tfunc,'pwl') == 1
            warning('Performing temporal adjustment with more than one tree can lead to oscillatory solution');
        end
    end
    if verbose == 1
        fprintf(1,'In %s, rank defiency persists! Consider reducing the number of tbreaks.\n',mfilename);
    end
end

%% Initialize Least Squares
if verbose == 1
    fprintf(1,'Begin least squares adjustment...\n');
end

%    ls_with_cov assumes that the covariance matrix of B is known only up to a
%    scale factor.  MSE is an estimate of that unknown scale factor, and
%    ls_with_cov scales the outputs S and STDX appropriately.  However, if V is
%    known to be exactly the covariance matrix of B, then that scaling is
%    unnecessary.  To get the appropriate estimates in this case, you should
%    rescale S and STDX by 1/MSE and sqrt(1/MSE), respectively.

if ndat == 1
    pest(1) = 0;  % arbitrarily make first epoch the origin
    pest(2) = data(1);
    psig(1) = 0;
    psig(2) = dsig(1);
    mse = 1;
else %continue with inversion
    A = G;
    B = colvec(data);
    
    %% Build covariance matrix
    %use incidence_to_cov to find the covariance matrix from edge
    %Laplacian; reference Biggs, et al. 2007 and Merris (1994)
    
    V = incidence_to_cov(Q, dsig);
    
    %     figure
    %     spy(V);
    %     title('data covariance matrix V');
    
    %% Perform the inversion using classic linear algebra solution to least-squares with pseudo-inverse
    
    if strcmp('ber_tikh', tfunc) == 1 && ncols == 1 % Case of Berardino et al. (2002) parameterization with Tikhonov regularization
        
        order = metaparams(1); % order for Tikhonov regularization
        beta_range = metaparams(2):metaparams(4):metaparams(3); % range for regularization parameter
        L2 = tikh_l_matrix(mparams,order);
        
        % if weighting the design matrix G use lsqr_tikh...
        [ pest, psig, mse, Vx, index] = lsqr_tikh( G, V, B, L2, beta_range);
        
    else   % All other parameterizations
        
        % Ginv is not specified, so we have to calculate it
        if dofast == 0
            Ginv = [];
        end
        [pest, psig, mse, Vx, Ginv] = ls_with_cov(A, B, V, Ginv);
        
        
        % display results
        if verbose == 1
            fprintf(1,'i,pest(i), psig(i)\n');
            for i=1:mparams
                fprintf (1,'%5d %#12.4e %#12.4e\n',i,pest(i),psig(i));
            end
        end
        
        
        %% Interpolate results to find displacement per epoch
        
        % modeled values for each pair
        dmod = G * pest;
        
        %dmod = interp1(tbreaks,pest,tm,'linear');
        dmod = dmod(1:ndat);
        
        % residuals
        res = data(1:ndat) - dmod;
        
        % consider only observations (not constraints)
        dsig = dsig(1:ndat);
        res_n = res./dsig; % residuals normalized by uncertainties
        
        % integrate over time
        if strcmp(tfunc, 'pwl') == 1 && ntrees == 1
            % for this system, the number of parameters equals the number of intervals + 1,
            % so we have no estimate of uncertainty
            tfit = tu;
            pfit = zeros(size(tfit));
            dm = 0;
            for i=2:numel(pest)
                dm = dm + pest(i) - pest(i-1);
                pfit(i) = dm; % estimated modeled displacement
            end
            sfitl = pfit;
            sfitu = pfit;
        else
            % all other functions use interpolation to convert from temporal
            % adjustment
            nfit = 1000;
            tfit1 = colvec(linspace(nanmin(tu),nanmax(tu),nfit));
            tfit = unique(sort(cat(1,tfit1,tu)));
            pfit = zeros(size(tfit));
            sfitl = zeros(size(tfit));
            sfitu = zeros(size(tfit));
            switch tfunc
                case {'nseg0','forge'} % need t2nd for these time functions
                    t2nd = nanmin(tu);
                otherwise
                    t2nd = nan;
            end
            for i=1:numel(tfit)
                %tf = time_function(tfunc, tfit(i),   tbreaks, metaparams); % find time value for given parameterization
                tf = time_function(tfunc, tfit(i),   tbreaks, metaparams,t2nd); % find time value for given parameterization
                dm=0;
                dsl = 0;
                dsu = 0;
                for j=1:mparams
                    dm  = dm   + tf(j) * pest(j); % interpolate to find displacement per epoch
                    dsl = dsl + tf(j) * (pest(j)-psig(j));
                    dsu = dsu + tf(j) * (pest(j)+psig(j));
                end
                pfit(i) = dm; % estimated model displacement
                sfitl(i) = dsl; % lower 1-sigma uncerainty bounds
                sfitu(i) = dsu; % upper 1-sigma uncertainty bounds
            end % loop over epochs
        end
    end
    
    %% Print results
    if verbose == 1
        fprintf(1,   'index  dateM dateS observed modeled residual res/sigma\n');
        for i = 1:ndat
            fprintf (1,'%5d %#12.4f %#12.4f %#12.4f %#12.4f %#12.4f %#12.4f %#10.4f\n',i,tm(i),ts(i),data(i),dmod(i),res(i),data_sigma(i),res(i)./data_sigma(i));
        end
    end
    
    % calculate chi-squared ourselves
    sswr = sum((res./dsig).^2)/numel(res);
    var=sum((res).^2)/(numel(res)-mparams);
    chisquare=sum((res./dsig).^2)/(numel(res));
    wvar=sum((res./dsig).^2)/(numel(res)-mparams);
    
    if verbose == 1
        fprintf(1,'MSE,       mean standard error (or factor used to scale variance of estimated parameters) is %20.4g\n',mse);
        fprintf(1,'SSWR,      sum of squared weighted residuals                                              is %20.4g\n',sswr);
        fprintf(1,'sqrt(MSE), square root of MSE                                                             is %20.4g\n',sqrt(mse));
        fprintf(1,'sqrt(SSWR),square root SSWR                                                               is %20.4g\n',sqrt(sswr));
        fprintf(1,'var, Sample variance in units of data                                                     is %20.4g\n',var);
        fprintf(1,'sqrt(var), Sample standard deviation in units of data                                     is %20.4g\n',sqrt(var));
        fprintf(1,'chisquare, Chisquared                                                                     is %20.4g\n',chisquare);
        fprintf(1,'wvar, estimate of normalized sample variance                                              is %20.4g\n',wvar);
        if strcmp('ber_tikh', tfunc) == 1
            fprintf(1,'index,     index of Tikhonov regularization parameter                                     is %20.4g\n',index);
            fprintf(1,'beta,      Tikhonov regularization parameter                                              is %20.4g\n',beta_range(index));
        end
    end
end
return
end



