function [pest, psig, mse, dmod, tfit, pfit, sfitl, sfitu, rd, V, G, data, sswr, Vx, var] = temporal_adjustment(data,data_sigma,tm,ts,tbreaks,tfunc,metaparams)
%function [pest, psig, mse, dmod, tfit, pfit, sfitl, sfitu, rd] = temporal_adjustment7(data,data_sigma,tm,ts,tbreaks,tfunc,metaparams)
% given time tags and data values, perform temporal adjustment
% use step functions
% by epochs in tbreaks
%%dmod = modeled value of differential quantity (from 1st epoch (tm) to 2nd epoch (ts))
%    intmod = modeled value of integrated quantity since beginning of time
%    (minimum of [tm, ts]
% tfunc == time function
%          'step'
%          'pwl'
%          'poly'
% Kurt Feigl University of Wisconsin-Madison
% 2011-OCT-08
% 2014-07-05 version 3 to handle Tabrez' test case
% 2014-07-15 version 4 to handle time functions including:
%    'nsegs'
% 2014-0729  version 5 to consider full data covariance matrix
% 20140824 Version 7 - return rank defiency
% 2014-09-09 Version 8b - full data covariance matrix using graoh theory- Elena Baluyut

fprintf(1,'%s begins ...\n',mfilename);


if exist('metaparams','var') == 0
    metaparams = nan;
end

iok=which('lscov');
if numel(iok) < 1
    error('Need lscov function')
end

% number of data
ndat = numel(data);
if numel(data_sigma) == ndat
    dsig = colvec(data_sigma);
else
    error('data_sigma has wrong dimensions');
end
dmod = nan(size(data));

% make column vectors
tm = colvec(tm); % master epoch
ts = colvec(ts); % slave epoch
tu = colvec(sort(unique([tm; ts]))); % unique set of epochs
me = numel(tu);  % number of epochs
fprintf(1,'Number of unique epochs me = %d\n',me);

%if  min(tbreaks) - min(tu) > 1.0/365.25
if  min(tbreaks) <  min(tu)
    error('min(tbreaks) before first observation');
end

if max(tbreaks) > max(tu)
    error('max(tbreaks) after last observation');
end



%% find trees
dispflag = 1;
[trees, DD] = findtrees(tm,ts,dispflag);
[ntrees,ndummy] = size(trees);
pairs = pairlist(DD);

%save the original incidence matrix, to be used later to create
%correlation and covariance matrices
DD0 = DD;
dsig0 = dsig;
data0 = data;


% sort time epochs
tbreaks = colvec(sort(unique(tbreaks)));
fprintf(1,'i,tbreaks(i)\n');
for i=1:numel(tbreaks)
    fprintf(1,'%d %12.5f\n',i,tbreaks(i));
end


% % build design matrix

if strcmp(tfunc,'pwl') == 1 && ntrees == 1
    % In this case, we can use the incidence matrix DD as the desgin
    % matrix G with an additional constraint of first epoch at 0
    % displacement
    [ndummy,mparams] = size(DD);
    G = DD;
    G(end+1, 1) = 1;
    data(end+1) = 0;
    dsig(end+1) = 1;
    
else
    %For all other cases, use temporal adjustment
    % get number of parameters
    disp('number of parameters')
    mparams = numel(time_function(tfunc, min([tm; tu]), tbreaks, metaparams));
    G = zeros(ndat,mparams);
    
    for i=1:ndat
        tfm = time_function(tfunc, tm(i), tbreaks, metaparams);
        tfs = time_function(tfunc, ts(i), tbreaks, metaparams);
        for j=1:mparams
            G(i,j) = tfs(j)-tfm(j);
        end
    end % loop over data
    
end

%Print Data Summary
fprintf(1,'i,data(i),dsig(i)\n');
for i = 1:numel(data)
    fprintf (1,'%5d %12.4g %12.4g \n',i,data(i),dsig(i));
end

disp 'Length of data vector, excluding constraints';ndat
disp 'Dimensons of design matrix, including constraints';nm = size(G)
disp 'Rank of G'; rank(G)
disp 'Rank defiency, including constraints'; rd = mparams - rank(G)
if rd > 0
    if ntrees > 1
        warning(sprintf('Number of trees (%d) is greater than 1.\n',ntrees));
        if strcmp(tfunc,'pwl') == 1
            warning('Performing temporal adjustment with more than one tree leads to incorrect oscillatory solution');
        end
    end
    figure;spy(G);hold on;
    xlabel('column');ylabel('row');
    title(sprintf('G matrix for %s has ndat = %d rows and mparams = %d columns'...
        ,strrep(mfilename,'_',' '),ndat,mparams));
    
    fprintf(1,'Consider reducing the number of tbreaks.\n');
    warning('Rank defiency persists!');
end

disp 'Begin least squares adjustment...'
%    LSCOV assumes that the covariance matrix of B is known only up to a
%     scale factor.  MSE is an estimate of that unknown scale factor, and
%     LSCOV scales the outputs S and STDX appropriately.  However, if V is
%     known to be exactly the covariance matrix of B, then that scaling is
%     unnecessary.  To get the appropriate estimates in this case, you should
%     rescale S and STDX by 1/MSE and sqrt(1/MSE), respectively.

if ndat < 2    pest(1) = 0;  % arbitrarily make first epoch the origin
    pest(2) = data(1);
    psig(1) = 0;
    psig(2) = dsig(1);
    mse = 1;
else
    A = G;
    B = colvec(data);
    
    %build covariance matrix
    %use incidence_to_cov to find the covariance matrix from edge Laplacian; reference Biggs, et al. 2007
    
    
    V = incidence_to_cov(DD, dsig, tfunc, ndat);
    
    figure
    spy(V);
    title('data covariance matrix V');
    
    
    % perform the inversion using classic linear algebra solution to
    % least-squares with pseudo-inverse
    
    [pest psig mse Vx] = ls_with_cov(A, B, V);
    
    fprintf(1,'i,pest(i), psig(i)\n');
    for i=1:mparams
        fprintf (1,'%5d %#12.4e %#12.4e\n',i,pest(i),psig(i));
    end
    
    
    % modeled values for each pair
    dmod = G * pest;
    
    %dmod = interp1(tbreaks,pest,tm,'linear');
    dmod = dmod(1:ndat);
    
    
    
    % residuals
    res = data(1:ndat) - dmod;
    
    
    % consider only observations (not constraints)
    dsig = dsig(1:ndat);
    res_n = res./dsig;
    
    % integrate over time
    if strcmp(tfunc, 'pwl') == 1 && ntrees == 1
        % for this system, the number of parameters equals the number of intervals + 1,
        % so we have no estimate of uncertainty
        tfit = tu;
        pfit = zeros(size(tfit));
        dm = 0;
        for i=2:numel(pest)
            dm = dm + pest(i) - pest(i-1);
            pfit(i) = dm;
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
        for i=1:numel(tfit)
            tf = time_function(tfunc, tfit(i),   tbreaks, metaparams);
            dm=0;
            dsl = 0;
            dsu = 0;
            for j=1:mparams
                dm  = dm   + tf(j) * pest(j);
                dsl = dsl + tf(j) * (pest(j)-psig(j));
                dsu = dsu + tf(j) * (pest(j)+psig(j));
            end
            pfit(i) = dm;
            sfitl(i) = dsl;
            sfitu(i) = dsu;
        end % loop over epochs
    end
end

fprintf(1,   'index  dateM dateS observed modeled residual res/sigma\n');
for i = 1:ndat
    fprintf (1,'%5d %#12.4f %#12.4f %#12.4g %#12.4g %#12.4g %#12.4g %#10.4f\n',i,tm(i),ts(i),data(i),dmod(i),res(i),data_sigma(i),res(i)./data_sigma(i));
end

% calculate chi-squared ourselves
sswr = sum((res./dsig).^2)/numel(res);
var=sum((res).^2)/(numel(res)-mparams);
chisquare=sum((res./dsig).^2)/(numel(res));
wvar=sum((res./dsig).^2)/(numel(res)-mparams);
fprintf(1,'MSE,       mean standard error (or factor used to scale variance of estimated parameters) is %20.4f\n',mse);
fprintf(1,'SSWR,      sum of squared weighted residuls                                               is %20.4f\n',sswr);
fprintf(1,'sqrt(MSE), square root of MSE                                                             is %20.4f\n',sqrt(mse));
fprintf(1,'sqrt(SSWR),square root SSWR                                                               is %20.4f\n',sqrt(sswr));
fprintf(1,'var, Sample variance in units of data                                                     is %20.4e\n',var);
fprintf(1,'sqrt(var), Sample standard deviation in units of data                                     is %20.4e\n',sqrt(var));
fprintf(1,'chisquare, Chisquared                                                                     is %20.4f\n',chisquare);
fprintf(1,'wvar, estimate of normalized sample variance                                              is %20.4f\n',wvar);
return



