function [pest, psig, mse, dmod, intmod, tfit, pfit, sfit,rfit] = adjust_steps3(data,data_sigma,tm,ts,imast,islav,tbreaks,DD)
%function [pest, psig, mse, dmod, dmodsig, tfit, pfit, sfit] = adjust_steps(data,data_sigma,tm,ts,imast,islav,tbreaks,DD)
% given time tags and data values, perform temporal adjustment
% use step functions
% by epochs in tbreaks
%%dmod = modeled value of differential quantity (from 1st epoch (tm) to 2nd epoch (ts))
%    intmod = modeled value of integrated quantity since beginning of time
%    (minimum of [tm, ts]
% Kurt Feigl University of Wisconsin-Madison
% 2011-OCT-08

fprintf(1,'%s begins ...\n',mfilename);

error(nargchk(8,8,nargin));
error(nargoutchk(9,9,nargout));

iok=which('lscov');
if numel(iok) < 1
    error('Need lscov function')
end

% number of data
ndat = length(data);
if numel(data_sigma) == ndat
    dsig = colvec(data_sigma);
else
    error('data_sigma has wrong dimensions');
end
dmod = nan(size(data));

% make column vectors
tm = colvec(tm); % master epoch
ts = colvec(ts); % slave epoch
tu = colvec(unique([tm; ts])); % unique set of epochs
me = numel(tu);  % number of epochs
[nrows,mcols] = size(DD);

% check number of epochs
if mcols ~= me
    error(sprintf('miscount in columns: found %d expected: %d\n',mcols,me));
end
% check number of pairs
if nrows ~= ndat
    error(sprintf('miscount in rows: found %d expected: %d\n',nrows,ndat));
end


% break points  
tbreaks = colvec(unique(tbreaks));


disp 'number of parameters equals number of break points';
mparam = numel(tbreaks) % number of break points

%G = zeros(ndat,me-1);
G = zeros(ndat,mparam);

for i=1:ndat
   fms = time_function4('step',tm(i),ts(i),tbreaks);
   for j=1:mparam
        G(i,j) = fms(j);
    end
end % loop over data

% length of intervals
dt = diff(tbreaks);

% fprintf(1,'i,j,G(i,j)\n');
% for i=1:ndat
%     for j=1:mparam
%         fprintf (1,'%5d %5d %12.4g\n',i,j,G(i,j));
%     end
% end

figure;spy(G);hold on;
xlabel('column');ylabel('row');
title(sprintf('G matrix for %s',strrep(mfilename,'_',' ')));

fprintf(1,'i,data(i),dsig(i)\n');
for i = 1:ndat
    fprintf (1,'%5d %12.4g %12.4g \n',i,data(i),dsig(i));
end

disp 'Length of data vector, including constraints';ndat = numel((data))
disp 'Dimensons of design matrix, including constraints';nm = size(G)
disp 'Rank defiency, including constraints'; rd = mparam - rank(G)
if rd > 0
    warning('Rank defiency persists!');
end

disp 'Begin least squares adjustment...'
%    LSCOV assumes that the covariance matrix of B is known only up to a
%     scale factor.  MSE is an estimate of that unknown scale factor, and
%     LSCOV scales the outputs S and STDX appropriately.  However, if V is
%     known to be exactly the covariance matrix of B, then that scaling is
%     unnecessary.  To get the appropriate estimates in this case, you should
%     rescale S and STDX by 1/MSE and sqrt(1/MSE), respectively.
 
if ndat < 2
    pest(1) = 0;  % arbitrarily make first the origin
    pest(2) = data(1);
    psig(1) = 0;
    psig(2) = dsig(1);
    mse = 1;
else
    A = G;
    B = colvec(data);
    W = colvec(1./dsig.^2);
    %W = eye(ndat); % unweighted least squares
    [pest,psig,mse]=lscov(A,B,W,'orth');
 
    % modeled values for each pair
    dmod = G * pest;
    %dmod = interp1(tbreaks,pest,tm,'linear');
    
    % integrate over time
    intmod = nan(size(data));
    %tfit = colvec(unique(sort([colvec(tu); colvec(tbreaks);colvec(tbreaks) - 1.0/365.25])));
    %tfit = colvec(linspace(floor(nanmin(tu)),ceil(nanmax(tu)),100));
    %tfit = colvec(unique(sort([colvec(tu); colvec(tbreaks)])));
    tfit = colvec(linspace(nanmin(tu),nanmax(tu),1000));
  

    pfit = zeros(size(tfit));
    sfit = zeros(size(tfit));
    rfit = nan(size(sfit));
    %ds = 0;
    for i=1:numel(tfit)
        dm = 0;
        ds2 = 0;
         for j = 1:mparam
            dt = tfit(i)-tbreaks(j);                         % time increment
            dm = dm + heaviside1(dt) * pest(j);              % accumulate modeled displacement
            ds2 = ds2 + (heaviside1(dt) * psig(j))^2;        % accumulate uncertainty in quadrature
        end
        pfit(i) = dm;
        sfit(i) = sqrt(ds2);        
    end % loop over data
end

fprintf(1,'i,tbreaks(i), pest(i), psig(i)\n');
for i=1:mparam
    fprintf (1,'%5d %#12.4f %#12.4e %#12.4e\n',i,tbreaks(i),pest(i),psig(i));
end

% tell us about the individual epochs
fprintf(1,   'index orbnum date     species yr Adjusted Sigma\n');
for i = 1:mparam
    fprintf (1,'%5d %#12.4f %#12.4g %#12.4g\n',i,tbreaks(i),pest(i),psig(i));
end

res = data-dmod;
%rfit=res;
ndat=numel(res)
nparam=4;

fprintf(1,   'index  dateM dateS observed modeled residual res/sigma\n');
for i = 1:ndat
    fprintf (1,'%5d %#12.4f %#12.4f %#12.4g %#12.4g %#12.4g %#12.4g %#10.4f\n',i,tm(i),ts(i),data(i),dmod(i),res(i),data_sigma(i),res(i)./data_sigma(i));
end

% calculate chi-squared ourselves
sswr = sum((res./dsig).^2)/numel(res);
var=sum((res).^2)/(numel(res)-nparam);
chisquare=sum((res./dsig).^2)/(numel(res));
wvar=sum((res./dsig).^2)/(numel(res)-nparam);
fprintf(1,'MSE,       mean standard error (or factor used to scale variance of estimated parameters) is %20.4f\n',mse);
fprintf(1,'SSWR,      sum of squared weighted residuls                                               is %20.4f\n',sswr);
fprintf(1,'sqrt(MSE), square root of MSE                                                             is %20.4f\n',sqrt(mse));
fprintf(1,'sqrt(SSWR),square root SSWR                                                               is %20.4f\n',sqrt(sswr));
fprintf(1,'var, Sample variance in units of data                                                     is %20.4e\n',var);
fprintf(1,'sqrt(var), Sample standard deviation in units of data                                     is %20.4e\n',sqrt(var));
fprintf(1,'chisquare, Chisquared                                                                     is %20.4f\n',chisquare);
fprintf(1,'wvar, estimate of normalized sample variance                                              is %20.4f\n',wvar);
return



