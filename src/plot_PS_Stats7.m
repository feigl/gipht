% plot statistics for Prolate Spheroid Model
% Kurt Feigl, 2011-06-13
% Summer, 2011-07-12 now plotting rates in fig 1
% Kurt 2011-09-22 Kurt, for NASA proposal
% 2011-10-05 Kurt use PST structure
% Summer 2011-11-22 now plots cases WITHOUT ensemble estimates
clear all; close all; fclose('all');
nf=0; % number of figures

%%% input: change file names by dataset  %%%
my_pnames={ 'time_fn_@_epoch_001_in_years____'...
    ,'time_fn_@_epoch_002_in_years____'...
    ,'YangPS_Easting_in_m_____________'...
    ,'YangPS_Northing_in_m____________'...
    ,'YangPS_Depth_in_m_______________'...
    ,'YangPS_Excess_Pressure_in_Pa____'...
    ,'YangPS_semimajor_axis_a__m______'...
    ,'YangPS_semiminor_axis_b_in_m____'...
    ,'YangPS_Azimuth_deg_CCW_from_N___'...
    ,'YangPS_Plunge_in_degrees________'};


% list of parameter files
% fnames = flist('results_T72_2/p*/PST.OUT');
% dataset = 'T72_2';
%fnames = flist('/data/okmok/gipht21/p_*x_20111026_*/PST.OUT');
%dataset = 'T344c_forUSGSREPORT';
%fnames = flist('/data/okmok/gipht21/results_summer/samp_gt_2000/*PST.OUT')
fnames = flist('/data/okmok/gipht21/results_summer/*T*/p*/PST.OUT') %

nfiles = numel(fnames)
fnames{end}

% add two co-eruptives
% fnames = strcat(fnames,{'/data/okmok/gipht22/KEEPx_20111203_112652/PST.OUT'}); % 1997
% fnames = strcat(fnames,{'/data/okmok/gipht22b/KEEPx_20111203_112740/PST.OUT'}); % 2008
% fnames{end+1} = '/data/okmok/gipht22/KEEPx_20111203_112652/PST.OUT'; % 1997
% fnames{end+1} = '/data/okmok/gipht22b/KEEPx_20111203_112740/PST.OUT'; % 2008
% 
% nfiles = numel(fnames)
% fnames{end}

dataset = 'InterPlus2Co';


my_pnames={ 'time_fn_@_epoch_001_in_years____'...
    ,'time_fn_@_epoch_002_in_years____'...
    ,'YangPS_Easting_in_m_____________'...
    ,'YangPS_Northing_in_m____________'...
    ,'YangPS_Depth_in_m_______________'...
    ,'YangPS_Excess_Pressure_in_Pa____'...
    ,'YangPS_semimajor_axis_a__m______'...
    ,'YangPS_semiminor_axis_b_in_m____'...
    ,'YangPS_Azimuth_deg_CCW_from_N___'...
    ,'YangPS_Plunge_in_degrees________'};

mparam = numel(my_pnames)

% divide by this number to convert from SI units (e.g., Pa)
% to user-friendly units (e.g., MPa)
%scalefactor = 1.0e6;
scalefactor = 1.0e3; % convert meters to km

pname1 = 'YangPS_Excess_Pressure_in_Pa____';
%pname1 = 'Pair_00001_Unwrap_Obs_Mod_mm____';

%titlestring = sprintf('Mogi Model %s\n',strrep(dataset,'_',''));
titlestring = sprintf('Yang Prolate Spheroid Model %s\n',strrep(dataset,'_',''));
%ylab = 'Pressure (MPa)';
ylab = 'Range Change at Center (mm)';

% retrieve the values of the parameters from the PST files
[P, pnames] = get_psts(fnames);
npairs = numel(P)


if npairs ~= nfiles
    error('miscount');
end
Pp1={P.p1};
Pps={P.sigma};

% field names are defined in this order
Yang.m  = nan(npairs,1);
Yang.s  = nan(npairs,1);
Yang.x  = nan(npairs,1);
Yang.y  = nan(npairs,1);
Yang.z  = nan(npairs,1);
Yang.p  = nan(npairs,1);
Yang.a  = nan(npairs,1);
Yang.b  = nan(npairs,1);
Yang.c  = nan(npairs,1);
Yang.d  = nan(npairs,1);
Siga.m  = nan(npairs,1);
Siga.s  = nan(npairs,1);
Siga.x  = nan(npairs,1);
Siga.y  = nan(npairs,1);
Siga.z  = nan(npairs,1);
Siga.p  = nan(npairs,1);
Siga.a  = nan(npairs,1);
Siga.b  = nan(npairs,1);
Siga.c  = nan(npairs,1);
Siga.d  = nan(npairs,1);

for i=1:npairs
    Yang.m(i)  = Pp1{i}(get_parameter_index(my_pnames{ 1},pnames)); % 'time_fn_@_epoch_001_in_years____'
    Yang.s(i)  = Pp1{i}(get_parameter_index(my_pnames{ 2},pnames)); % 'time_fn_@_epoch_002_in_years____'
    Yang.x(i)  = Pp1{i}(get_parameter_index(my_pnames{ 3},pnames)); % 'YangPS_Easting_in_m_____________'
    Yang.y(i)  = Pp1{i}(get_parameter_index(my_pnames{ 4},pnames)); % 'YangPS_Northing_in_m____________'
    Yang.z(i)  = Pp1{i}(get_parameter_index(my_pnames{ 5},pnames)); % 'YangPS_Depth_in_m_______________'
    Yang.p(i)  = Pp1{i}(get_parameter_index(my_pnames{ 6},pnames)); % 'YangPS_Excess_Pressure_in_Pa____'
    Yang.a(i)  = Pp1{i}(get_parameter_index(my_pnames{ 7},pnames)); % 'YangPS_semimajor_axis_a__m______'
    Yang.b(i)  = Pp1{i}(get_parameter_index(my_pnames{ 8},pnames)); % 'YangPS_semiminor_axis_b_in_m____'
    Yang.c(i)  = Pp1{i}(get_parameter_index(my_pnames{ 9},pnames)); % 'YangPS_Azimuth_deg_CCW_from_N___'
    Yang.d(i)  = Pp1{i}(get_parameter_index(my_pnames{10},pnames)); % 'YangPS_Plunge_in_degrees________'
    %     Siga.m(i)  = Pps{i}(get_parameter_index(my_pnames{ 1},pnames));
    %     Siga.s(i)  = Pps{i}(get_parameter_index(my_pnames{ 2},pnames));
    Siga.m(i)  = NaN;
    Siga.s(i)  = NaN;
    Siga.x(i)  = Pps{i}(get_parameter_index(my_pnames{ 3},pnames));
    Siga.y(i)  = Pps{i}(get_parameter_index(my_pnames{ 4},pnames));
    Siga.z(i)  = Pps{i}(get_parameter_index(my_pnames{ 5},pnames));
    Siga.p(i)  = Pps{i}(get_parameter_index(my_pnames{ 6},pnames));
    Siga.a(i)  = Pps{i}(get_parameter_index(my_pnames{ 7},pnames));
    Siga.b(i)  = Pps{i}(get_parameter_index(my_pnames{ 8},pnames));
    Siga.c(i)  = Pps{i}(get_parameter_index(my_pnames{ 9},pnames));
    Siga.d(i)  = Pps{i}(get_parameter_index(my_pnames{10},pnames));
end

field_names=fieldnames(Yang)

% derived quantities
Ense.m = min([Yang.m Yang.s]);
Ense.s = max([Yang.m Yang.s]);
Sige.m = nan(size(Ense.m));
Sige.s = nan(size(Ense.s));
tmide = mean((Yang.m+Yang.s)/2.);
dte   = abs(Yang.m-Yang.s);
Yang.tmid = colvec((Yang.m + Yang.s)/2.); % middle of time interval
Yang.dt   = colvec(abs(Yang.s - Yang.m)); % duration of time interval
Yang.dataset = dataset

Sstd.name = 'standard deviation';
Wavg.name = 'weighted average';
Wstd.name = 'weighted std';
Sfac.name = 'scale factor for std';
Sigb.name = 'rescaled sigma';
Savg.name = 'unweighted mean';

for i = 1:mparam
    xvals = getfield(Yang,field_names{i});
    sigas = getfield(Siga,field_names{i});
    % estimate unweighted mean
    xmen = nanmean(xvals)
    Savg = setfield(Savg,field_names{i},xmen);
    % estimate sample standard deviations
    sstd = nanstd(xvals);
    Sstd = setfield(Sstd,field_names{i},sstd);
    % estimate weighted means, assuming each parameter is independent
    [wavg,wstd,sfac,sigb] = wmean(xvals,sigas);
    Wavg = setfield(Wavg,field_names{i},wavg);
    Wstd = setfield(Wstd,field_names{i},wstd);
    Sfac = setfield(Sfac,field_names{i},sfac);
    Sigb = setfield(Sigb,field_names{i},sigb);
end
% consider adding dataset to this filename
save(sprintf('%s_%s.mat',mfilename,dataset),'pnames','P','Yang','Ense','Siga','Sige','Sfac','Sigb'); %'pnames','E','P'

return


% prune out insignificant rates
% sigrat = abs(pr) ./ rs;
% iok = find(sigrat > 2.0 & isfinite(sigrat) == 1);
% prune out NaN
% iok = find(isfinite(rs) == 1);
% tm = colvec(tm(iok)); % master epoch in years
% Yang.s = colvec(Yang.s(iok)); % slave epoch in years
% pr = colvec(pr(iok)); % annual rate of change in pressure (pa/yr)
% rs = colvec(rs(iok)); % uncerainty (pa/yr)
% ibad = find(isnan(rs));
% rs(ibad)=1;

ndat = numel(Yang.m)

% list and plot estimated values
fidout = fopen(sprintf('%s_%s.out',mfilename,dataset),'w');
for fid = [1 fidout]
    for j = 1:mparam
        % list estimated values
        nf=nf+1;h(nf)=figure;
        subplot(2,1,1);axis ij;
        axis([1,50,-2,ndat + 4]);
        hold off;
        axis off;
 
        fmt='%5s %12.4e %12.4e\n';
        text(1,0,sprintf('%s %s',strrep(field_names{j},'_',' '),strrep(my_pnames{j},'_',' '))...
            ,'HorizontalAlignment','Left'...
            ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
        %ense = getfield(Ense,field_names{j}); % comment out when no ensemble
        % sige = getfield(Sige,field_names{j}); % comment out when no ensemble
        vals = getfield(Yang,field_names{j});
        siga = getfield(Siga,field_names{j});
        wavg = getfield(Wavg,field_names{j});
        wstd = getfield(Wstd,field_names{j});
        sfac = getfield(Sfac,field_names{j});
        savg = getfield(Savg,field_names{j});
        sstd = getfield(Sstd,field_names{j});
        for k = 1:ndat
            fprintf(fid,fmt,sprintf('%5d',k),vals(k),siga(k));
            text(1,k,sprintf(fmt,sprintf('%5d',k),vals(k),siga(k))...
                ,'HorizontalAlignment','Left'...
                ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
        end
        %fprintf(fid,fmt,sprintf('%5s','Ense'),ense,sige); % comment out when no ensemble
        fprintf(fid,fmt,sprintf('%5s','Savg'),savg,sstd);
        fprintf(fid,fmt,sprintf('%5s','Wavg'),wavg,wstd);
        fprintf(fid,fmt,sprintf('%5s','Sfac'),nan,sfac);
        %k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Ense'),ense,sige),'HorizontalAlignment','Left'...
            %    ,'FontName','Courier','Fontsize',7,'FontWeight','normal');  % comment out when no ensemble
        k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Savg'),savg,sstd),'HorizontalAlignment','Left'...
                ,'FontName','Courier','Fontsize',7,'FontWeight','normal'); 
        k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Wavg'),wavg,wstd),'HorizontalAlignment','Left'...
                ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
        k=k+1;text(1,k,sprintf(fmt,sprintf('%5s','Sfac'),nan, sfac),'HorizontalAlignment','Left'...
                ,'FontName','Courier','Fontsize',7,'FontWeight','normal');
        h=title(sprintf('%s %s'...
            ,strrep(dataset,'_',' ')...
            ,strrep(field_names{j},'_',' ')));
           set(h,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');

        % plot estimated values      
        subplot(2,1,2);axis xy;hold on;

        errorbar_plus2(Yang.tmid,vals,        Yang.dt/2.0, siga,'b+', 5); % ,1 individual estimates in blue
        errorbar_plus2(tmide,    wavg,dte/2.0,wstd,             'gs',10); % ,3 weighted mean in green
       % errorbar_plus2(tmide,    ense,        dte/2.0,     sige,'r*',10);
       % % ,3 Ensemble in red % comment out when no ensemble
 
        fixlabels('year','%.1f',strrep(my_pnames{j},'_',' '),'%.2g');
        
        printpdf(sprintf('%s_%s_param%02d.pdf',mfilename,dataset,j));
 
      end
end
fid
fclose(fidout);

% Prune data set
% plot outliers with error bars
% iset = find(abs(Yang.x-Ense.x) > 3.0*Sigb.x);
% iset = union(iset,find(abs(Yang.y-Ense.y) > 3.0*Sigb.y));
% iset = union(iset,find(abs(Yang.z-Ense.z) > 3.0*Sigb.z));
% iset = union(iset,find(abs(Yang.p-Ense.p) > 3.0*Sigb.p));

% Select all
iset = 1:npairs;


%%%%%%%%%%%%%%%%%%%%
% Plot all parameters pairwise
nf=nf+1;h(nf)=figure;

%  GPLOTMATRIX  Scatter plot matrix with grouping variable.
%     GPLOTMATRIX(X,Y,G) creates a matrix of scatter plots of the columns of
%     X against the columns of Y, grouped by G.  If X is P-by-M and Y is
%     P-by-N, GPLOTMATRIX will produce a N-by-M matrix of axes.  If you omit
%     Y or specify it as [], the function graphs X vs. X.  G is a grouping
%     variable that determines the marker and color assigned to each point in
%     each matrix, and it can be a categorical variable, vector, string
%     matrix, or cell array of strings.  Alternatively G can be a cell array
%     of grouping variables (such as {G1 G2 G3}) to group the values in X by
%     each unique combination of grouping variable values.

XVALS = zeros(npairs,mparam);
YVALS = zeros(npairs,mparam);
for i=1:mparam % rows
    vals = colvec(getfield(Yang,field_names{i}));
    vals = vals-nanmean(vals);
    YVALS(:,i) = vals;
    XVALS(:,i) = vals;
end
disp 'XVALS';size(XVALS)
disp 'YVALS';size(YVALS)

% GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ,DOLEG,DISPOPT,XNAM,YNAM) specifies XNAM
%     and YNAM as the names of the X and Y variables.  Each must be a
%     character array or cell array of strings of the appropriate dimension.
G = 1:npairs;
doleg = 'off';
%dispopt = 'variable';
dispopt = 'on';
gplotmatrix(XVALS,YVALS,G...
    ,'krb','o+^v',3,doleg,dispopt...
    ,field_names,field_names);
h=title(strrep(dataset,'_',' '));set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
printpdf(sprintf('%s_%s_%02d.pdf',mfilename,dataset,nf));

%%%%%%%%%%%%%%%%%%%%
% Plot easting and northing in map view
nf=nf+1;h(nf)=figure;
subplot(2,1,1);
axis ij;hold on;axis equal;axis([-Inf +Inf -Inf Inf]);axis tight;
set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
% errorbar_plus2(Yang.x/scalefactor,Yang.y/scalefactor...
%               ,Siga.x/scalefactor,Siga.y/scalefactor...
%               ,'bx',1);
plot(Yang.x/scalefactor,Yang.y/scalefactor...
    ,'b+','MarkerSize',5);

% make prolate spheroid by repeating b axis
% [Ense.xs,Ense.ys,Ense.zs] =
% get_ellipsoid(Ense.x,Ense.y,Ense.z,Ense.a,Ense.b,Ense.b,Ense.c,Ense.d); % comment out when no ensemble

% draw surface projection of spheroid
%islice = find(abs(Ense.zs-Ense.z)<5.); % comment out when no ensemble
%plot(Ense.xs(islice)/scalefactor,Ense.ys(islice)/scalefactor,'r.','MarkerSize',1); % comment out when no ensemble

if numel(iset) > 0
    errorbar_plus2(Yang.x(iset)/scalefactor,Yang.y(iset)/scalefactor...
        ,Sigb.x(iset)/scalefactor,Sigb.y(iset)/scalefactor...
        ,'bx',5);
end
errorbar_plus2(Wavg.x/scalefactor,Wavg.y/scalefactor...
    ,Wstd.x/scalefactor,Wstd.y/scalefactor...
    ,'gs',10,3);
errorbar_plus2(Savg.x/scalefactor,Savg.y/scalefactor...
    ,Sstd.x/scalefactor,Sstd.y/scalefactor...
    ,'ko',10,3);
% errorbar_plus2(Ense.x/scalefactor,Ense.y/scalefactor...
  %  ,Sige.x/scalefactor,Sige.y/scalefactor...
  %  ,'r*',10,3); % comment out when no ensemble

%fixlabels(strrep('Yang1_Easting_in_km','_',' '),'%.0f',strrep('Yang1_Northing_in_km','_',' '),'%.1f');
fixlabels(strrep('Yang1_Easting_in_km','_',' '),'',strrep('Yang1_Northing_in_km','_',' '),'');
h=title(strrep(dataset,'_',' '));set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%
% Draw Cross Section along Y = Ycenter
subplot(2,1,2);
%nf=nf+1;h(nf)=figure;hold on;
axis ij;hold on;axis equal;axis([-Inf +Inf 0 Inf]);axis tight;

set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
plot(Yang.x/scalefactor,Yang.z/scalefactor...
    ,'b+','MarkerSize',5);

%islice = find(abs(Ense.ys-Ense.y)<5.); % comment out when no ensemble
%plot(Ense.xs(islice)/scalefactor,Ense.zs(islice)/scalefactor,'r.','MarkerSize',1);% comment out when no ensemble
if numel(iset) > 0
    errorbar_plus2(Yang.x(iset)/scalefactor,Yang.z(iset)/scalefactor...
        ,Sigb.x(iset)/scalefactor,Sigb.z(iset)/scalefactor...
        ,'bx',5);
end
errorbar_plus2(Wavg.x/scalefactor,Wavg.z/scalefactor...
    ,Wstd.x/scalefactor,Wstd.z/scalefactor...
    ,'gs',10,3);
errorbar_plus2(Savg.x/scalefactor,Savg.z/scalefactor...
    ,Sstd.x/scalefactor,Sstd.z/scalefactor...
    ,'ko',10,3);
%errorbar_plus2(Ense.x/scalefactor,Ense.z/scalefactor...
  %  ,Sige.x/scalefactor,Sige.z/scalefactor...
  %  ,'r*',10,3); % comment out when no ensemble
%fixlabels(strrep('Yang1_Easting_in_km','_',' '),'%.0f',strrep('Yang1_Depth_in_km','_',' '),'%.0f');
fixlabels(strrep('Yang1_Easting_in_km','_',' '),'',strrep('Yang1_Depth_in_km','_',' '),'');
%title(strrep(dataset,'_',' '));
%legend('individual pairs','outliers +/- scaled \sigma','weighted mean of ind.','mean','ensemble')
printpdf(sprintf('%s_%s_%02d.pdf',mfilename,dataset,nf));


%%%%%%%%%%%%%%%%%%%%
% Draw Pressure as a function of time
nf=nf+1;h(nf)=figure;hold on;
set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold');
plot(Yang.tmid,Yang.p/1.e6...
    ,'b+','MarkerSize',5);
if numel(iset) > 0
    errorbar_plus2(Yang.tmid(iset),Yang.p(iset)/1.e6...
        ,Yang.dt(iset)/2.0,Sigb.p(iset)/1.e6...
        ,'bx',5);
end
errorbar_plus2(tmide,  Wavg.p/1.e6 ...
    ,dte/2.0,Wstd.p/1.e6 ...
    ,'gs',10,3);
errorbar_plus2(tmide,  Savg.p/1.e6 ...
    ,dte/2.0,Sstd.p/1.e6 ...
    ,'ko',10,3);
%errorbar_plus2(tmide,  Ense.p/1.e6 ...
 %   ,dte/2.0,Sige.p/1.e6 ...
  %  ,'r*',10,3); % comment out when no ensemble
%fixlabels(strrep('Year','_',' '),'%.1f',strrep('Yearly_rate_of_Change_Yang1_Excess_Pressure_in_MPa_per_year','_',' '),'%.0f');
fixlabels(strrep('Year','_',' '),'%.1f',strrep('Yearly_rate_of_Change_Yang1_Excess_Pressure_in_MPa_per_year','_',' '),'');
%title(strrep(dataset,'_',' '));
%legend('individual pairs','outliers +/- scaled \sigma','weighted mean of ind.','mean','ensemble')
printpdf(sprintf('%s_%s_%02d.pdf',mfilename,dataset,nf));

return






