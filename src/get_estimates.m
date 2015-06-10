function matfilename = get_estimates(varargin)
%function matfilename = get_estimates(fnames,dataset,ename,mindt,minsigratio)
% extract values and uncertainties of estimated parameters from GIPHT solutions in PST.OUT files
% Kurt Feigl, 2014-02-03

fnames  = varargin{1};
dataset = varargin{2};

% full file name for PST file from Ensemble solution
if nargin >= 3
    ename  = varargin{3};
else
    ename = '';
end

% minimum time interval in years
if nargin >= 4
    mindt = varargin{4};
else
    mindt = 0.;
end

% mininimum significance ratio
if nargin >= 5
    minsigratio  = varargin{5};
else
    minsigratio  = 0.;
end

% initialize
nf=0; % number of figures
nfiles = numel(fnames); % number of files

% retrieve the values of the parameters from the PST files
[P, pnames] = get_psts(fnames);
npairs = numel(P)

if npairs ~= nfiles
    error('miscount');
end

% count parameters
mparams0 = numel(P(1).p1)
for i=2:npairs
    mparams = numel(P(i).p1);
    if mparams ~= mparams0
        warning('error in pair %i\n',i);
    end
end

mparams

% master and slave epochs
itmast = get_parameter_index('time_fn_@_epoch_001_in_years____',pnames);
itslav = get_parameter_index('time_fn_@_epoch_002_in_years____',pnames);
tmast = nan(npairs,1);
tslav = nan(npairs,1);
for i=1:npairs
    tmast(i) = P(i).p1(itmast);
    tslav(i) = P(i).p1(itslav);
    % should compress list of short pairs here
    %if abs(tmast(i)-tslav(i)) >= mindt
end

% count adjusted parameters with a finite uncertainty
ifree=[]; % pointer to parameter lists
for j=1:mparams
    kount = 0;
    for i=1:npairs
        pa = P(i).p1(j)-P(i).p0(j); % adjustment
        ps = P(i).sigma(j);         % uncertainty
        if isfinite(ps)==1 && abs(pa/ps) > minsigratio
            kount = kount+1;
        end
    end
    if kount == npairs
        ifree(end+1) = j;
    end
end
mfree = numel(ifree)


% print some values
for j=1:numel(ifree)
    k = ifree(j);
    for i=1:npairs
        fprintf(1,'%32s %9.4f %9.4f %20.4e %20.4e\n',pnames{k},tmast(i),tslav(i),P(i).p1(k),P(i).sigma(k));
    end
end

% ;prune values and calculate statistics
pvals = nan(npairs,mfree);
sigas = nan(npairs,mfree);
sigbs = nan(npairs,mfree);
xmen  = nan(mfree,1);
sstd  = nan(mfree,1);
wstd  = nan(mfree,1);
wavg  = nan(mfree,1);
sfac  = nan(mfree,1);
ssem  = nan(mfree,1);
kk=0;
pnames_free = cell(mfree,32);

for j=1:mfree
    jp = ifree(j); % pointer into list of parameters
    kk = kk+1;     % pointer into compressed list of parameters
    %pnames_free{kk} = strrep(pnames{jp},'@','_');
    pnames_free{kk} = pnames{jp};
    for i=1:npairs
        pvals(i,kk) = P(i).p1(jp);
        sigas(i,kk) = P(i).sigma(jp);
    end
    
    % estimate unweighted mean
    xmen(kk) = nanmean(pvals(:,kk));
    
    % estimate sample standard deviation
    sstd(kk) = nanstd(pvals(:,kk));

    % standard error of mean
    ssem(kk) = nanstd(pvals(:,kk))/sqrt(numel(find(isfinite(pvals(:,kk)))));

    % estimate weighted means, assuming each parameter is independent
    [wavg(kk),wstd(kk),sfac(kk),sigbs(:,kk)] = wmean(pvals(:,kk),sigas(:,kk));
    
    fprintf(1,'%32s %9.4f %9.4f %20.4e %20.4e\n',pnames{k},tmast(i),tslav(i),xmen(kk),sstd(kk));
end

% get Ensemble solution
[E, Enames] = get_psts(ename);
itense = get_parameter_index('time_fn_@_epoch',Enames);
tens = E.p1(itense);
tmine = nanmin(tens);
tmaxe = nanmax(tens);


% derived quantities
tmin = min([tmast tslav]);
tmax = max([tmast tslav]);
tmid = colvec((tmast + tslav)/2.); % middle of time interval
tdel  = colvec(abs(tslav - tmast)); % duration of time interval

% Pest.name = 'final estimate';
% Sstd.name = 'standard deviation';
% Wavg.name = 'weighted average';
% Wstd.name = 'weighted std';
% Sfac.name = 'scale factor for std';
% Sigb.name = 'rescaled sigma';
% Savg.name = 'unweighted mean';


% save files for later use
matfilename = sprintf('%s_%s.mat',mfilename,dataset);
save(matfilename,'pnames_free','tmast','tslav','pvals','sigas','xmen','sstd','wavg','wstd');



% prune out insignificant rates
% sigrat = abs(pr) ./ rs;
% iok = find(sigrat > 2.0 & isfinite(sigrat) == 1);
% prune out NaN
% iok = find(isfinite(rs) == 1);
% tm = colvec(tm(iok)); % master epoch in years
% tslav = colvec(tslav(iok)); % slave epoch in years
% pr = colvec(pr(iok)); % annual rate of change in pressure (pa/yr)
% rs = colvec(rs(iok)); % uncerainty (pa/yr)
% ibad = find(isnan(rs));
% rs(ibad)=1;

ndat = npairs;

% list and plot estimated values
% fidout = fopen(sprintf('%s_%s.out',mfilename,dataset),'w');
% for fid = [1 fidout]
fid = 1;
for j = 1:mfree
%    for j=5:5
    % list estimated values
    
    % find parameter in list for Ensemble
    ipe = get_parameter_index(pnames_free{j},Enames);
    %fprintf(1,'%s\n%s\n',pnames_free{j},Enames{ipe});
    if ipe > 0
        ense = E.p1(ipe);       % Ensemble estimate
        sige = E.sigma(ipe);    % Ensemble uncertainty
    else
        ense = NaN;
        sige = NaN;
    end

    nf=nf+1;h(nf)=figure;
    subplot(2,1,1);axis ij;
    axis([1,50,-2,ndat + 4]);
    hold off;
    axis off;
    
    fmt  ='%5s %12.4e %12.4e\n';
    fmt2 ='%5s %12.4e %12.4e %9.4f %9.4f\n';
    titlestr = strrep(pnames_free{j},'_',' ');
    fprintf(1,'%s\n',titlestr);
    text(1,0,titlestr...
        ,'HorizontalAlignment','Left'...
        ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
    for k = 1:npairs
        fprintf(fid,fmt2,sprintf('%5d',k),pvals(k,j),sigas(k,j),tmast(k),tslav(k));
        text(1,k,sprintf(fmt2,sprintf('%5d',k),pvals(k,j),sigas(k,j) ...
            ,tmast(k),tslav(k)) ...
            ,'HorizontalAlignment','Left'...
            ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
    end
    
    fprintf(fid,fmt,sprintf('%5s','Wavg'),wavg(j),wstd(j));
    fprintf(fid,fmt,sprintf('%5s','Sfac'),nan,sfac(j));
    %         k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Savg'),savg,sstd),'HorizontalAlignment','Left'...
    %                 ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
    k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Wavg'),wavg(j),wstd(j)),'HorizontalAlignment','Left'...
        ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
    k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Sfac'),nan, sfac(j)),'HorizontalAlignment','Left'...
        ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
    fprintf(fid,fmt,sprintf('%5s','Ense'),ense,sige); % comment out when no ensemble
    %        fprintf(fid,fmt,sprintf('%5s','Savg'),savg,sstd);
    
    if isfinite(ense) == 1
        k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Ense'),ense,sige),'HorizontalAlignment','Left'...
            ,'FontName','Courier','Fontsize',7,'FontWeight','normal');  % comment out when no ensemble
    end
    
    h=title(sprintf('%s %s'...
        ,strrep(dataset,'_',' ')...
        ,strrep(pnames_free{j},'_',' ')));
    set(h,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
    
    % plot estimated values
    subplot(2,1,2);axis xy;hold on;
    h=title(strrep(pnames_free{j},'_',' '));
    set(h,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
    
    % draw the symbols
    errorbar_plus2(tmid,                     pvals(:,j),                  tdel/2.0, sigas(:,j),          'b+',5);  % individual estimates in blue
    errorbar_plus2((max(tmax)+min(tmin))/2.,    wavg(j),  (max(tmax)-min(tmin))/2., wstd(j),             'gs',10); % weighted mean in green
    errorbar_plus2((max(tmax)+min(tmin))/2.,    xmen(j),  (max(tmax)-min(tmin))/2., ssem(j),             'ko',10); % unweighted mean and its standard error in black
    errorbar_plus2((tmaxe+tmine)/2.,               ense,  (tmaxe-tmine)/2.0,        sige,                'r*',10); % Ensemble in red 
    
    
    %fixlabels('year','%.1f',strrep(my_pnames{j},'_',' '),'%.2g');
    
    printpdf(sprintf('%s_%s_param%02d.pdf',mfilename,dataset,j));
    
end
% fid
% fclose(fidout);

return
end

