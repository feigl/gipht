%% script for determining the best InSAR pairs for analysis
% for use with pairs in a single track
% finds MST based on maximizing coherence 
% need file with the following columns:
% master slave bperp master_Dopplercentroid slave_Dopplercentroid
% note: doppler centroid can be found in GMTSAR .PRM files (param fd1)
% Elena C. Reinisch, 20190320

% initialize
clear all; close all;
addpath(genpath('~/PoroTomo'))
addpath(genpath('~/gipht'))
addpath(genpath('~/MATLAB'))
% addpath(genpath('~/Desktop/archGraphTreeTA'))

% set satellite dependent parameters
b_critical = 1100; %m; critical baseline set for ENVI and ERS (https://pdfs.semanticscholar.org/7791/7825df294c2b14d4ed0ce846de1b616932cd.pdf) (https://earth.esa.int/pub/ESA_DOC/envisat_val_1202/proceedings/ASAR/06_holzner.pdf)
weight_b_crit = 1/2; % restrict to be within 1/4 of critical baseline based on topography of region (e.g. ftp://topex.ucsd.edu/pub/sandwell/Venus.../critical_baseline/InSAR_VENUS.doc)
c = 0; % minimum bias correlation value
b_temporal = 80; % days 
azimuthal_bandwidth = 1652; %1680; % 1652 Hz for ENVI (from .PRM)

% load data
datfile = 'ERS_T27_sanem_bp_dopp_new.txt'
VALS = readtable(datfile);

[ndat ndum] = size(VALS);
ind = ndat;
orb1 = table2array(VALS(1:ind,1));
orb2 = table2array(VALS(1:ind,2));
bperp = abs(table2array(VALS(1:ind,3))); 
mast_dop = abs(table2array(VALS(1:ind,4)));
slav_dop = abs(table2array(VALS(1:ind,5)));

tm = orb1;
ts = orb2;
 
% check for twins
twin_ind = find(tm - ts == 0);
tm(twin_ind) = [];
ts(twin_ind) = [];
orb1(twin_ind) = []; 
orb2(twin_ind) = []; 
bperp(twin_ind) = [];  
mast_dop(twin_ind) = [];
slav_dop(twin_ind) = [];

% check for 0 bperp
bp_ind = find(bperp == 0);
bperp(bp_ind) = nanmean(bperp);  

% check for nan Dfreq
dp_ind = find(isnan(mast_dop) == 1);
mast_dop(dp_ind) = nanmean(mast_dop);
dp_ind = find(isnan(slav_dop) == 1);
slav_dop(dp_ind) = nanmean(slav_dop);


[ndat ndum] = size(tm);
ind = ndat;

sigd = ones(size(bperp));

tu = sort(unique([tm ts]));

% convert dates to datetimes
tu_datetime = cal2datetime(tu);
tm_datetime = cal2datetime(tm);
ts_datetime = cal2datetime(ts);

% get indices of master and slave epochs in tu
tm_ind = zeros(size(tm));
ts_ind = zeros(size(ts));
for i = 1:numel(tm)
    tm_ind(i) = find(tu == tm(i));
    ts_ind(i) = find(tu == ts(i));
end

% find Q
[trees, DD, tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(tm,ts,1);

format long

% plot bp
Qdiff2 = adjustbp(tepochs,DD,bperp, trees,iuniqorbs, uniqdates);
cal_date = [month(cal2datetime(tepochs)) day(cal2datetime(tepochs))];

% plot 
tepochs = dyear(year(cal2datetime(tepochs)), month(cal2datetime(tepochs)), day(cal2datetime(tepochs)));
figure;
set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
set(gcf,'DefaultTextInterpreter','tex');
plotbp_mst(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]','Bperp','m'));


% %% find minimum spanning tree using written function
% % > good
% i_weight = bperp;
% i_epochs_start = tm;
% i_epochs_end = ts;
% % find minmum spanning tree using Prim's algorithm and scripts from
% % GraphTreeTA
% [mdum, tm_ind] = ismember(tm, tu);
% [mdum, ts_ind] = ismember(ts, tu);
% [trees, DD_full tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(i_epochs_start,i_epochs_end,0, tm_ind, tm, ts_ind, ts);
% [i_start_msf, i_end_msf] = plotMSF( DD_full, trees, (i_weight), 1, 0);
% Best_pairs = [tu(i_start_msf), tu(i_end_msf)]
% % < good

%% Define Weights
% use bp and indices of epochs
corr_bp = (1 - abs(bperp)/(b_critical*weight_b_crit)).*max(heaviside(b_critical*weight_b_crit- abs(bperp)), 1e-2);
%temporal_range = exp(abs(caldays(between(ts_datetime, tm_datetime, 'd')))/b_temporal);
%temporal_normalized = (temporal_range - (min(temporal_range)+.001))./(max(temporal_range) - min(temporal_range));
%corr_t = temporal_normalized.*max(min(min(heaviside(month(tm_datetime) - 3.5) ,(1- heaviside(month(tm_datetime) - 11.5))), min(heaviside(month(ts_datetime) - 3.5) ,(1- heaviside(month(ts_datetime) - 11.5)))), 1e-10);
corr_t = exp((caldays(between(ts_datetime, tm_datetime, 'd')))/b_temporal).*max(min(min(heaviside(month(tm_datetime) - 3.5) ,(1- heaviside(month(tm_datetime) - 11.5))), min(heaviside(month(ts_datetime) - 3.5) ,(1- heaviside(month(ts_datetime) - 11.5)))), 1e-2);
corr_dop = (1 - abs(mast_dop - slav_dop)./azimuthal_bandwidth).*max(heaviside(azimuthal_bandwidth - abs(mast_dop - slav_dop)), 1e-2);
corr_sum = abs(corr_bp.*corr_t.*corr_dop);
% W = exp(-corr_sum*10);
W0 = (1 - normalize(corr_sum, 'range'))+1;

orb1_ind = [];
orb2_ind = [];
for i = 1:numel(orb1)
    orb1_ind(i,1) = find(tu == orb1(i));
    orb2_ind(i,1) = find(tu == orb2(i));
end

% add extra pair (connecting last epoch to itself) to make complete graph
% but weight so that it is never chosen
% W(end+1) = max(W) + 1;
% orb1_ind(end+1) = orb2_ind(end);
% orb2_ind(end+1) = orb2_ind(end);

%% find minimum spanning tree using MATLAB function
W = zeros(numel(tu));
for i = 1:numel(orb1_ind)
   W(orb1_ind(i), orb2_ind(i)) =  W0(i);
end

% build graph
% DG = sparse([orb1_ind], [orb2_ind], W);
DG = sparse(W);
UG = tril(DG+DG');

% find minimum spanning tree
t1 = tic;
[best_pairs, nodes] = graphminspantree(UG);
toc(t1)
% get MST pairs
% best_bp = (nonzeros(best_pairs));
% [mdum1, best_bp_ind] = ismember(best_bp, W);
% Best_pairs = [orb1(best_bp_ind), orb2(best_bp_ind)]
% Best_pairs = zeros(numel(best_bp), 2);
% for i = 1:numel(best_bp)
%     [bi, bj] = find(W == best_bp(i));
%     Best_pairs(i,:) = [orb1(bi), orb2(bj)];
% end

[ii, ij, best_weight] = find(best_pairs);
best_bp_ind = sort([ii, ij], 2);
Best_pairs = [tu(best_bp_ind)]

best_pairs_sort = zeros(size(best_pairs));
for i = 1:numel(best_weight)
    best_pairs_sort(best_bp_ind(i,1), best_bp_ind(i,2)) = best_weight(i);
end
best_pairs_sort = sparse(best_pairs_sort);

% view MST
bg = (biograph(best_pairs_sort,datestr(cal2datetime(tu), 'yyyy-mmm-dd'),'ShowArrows','on','ShowWeights','off'));
view(bg)

%dlmwrite('best_pairs_mst_ers27_bp.txt', Best_pairs, 'delimiter', '\t', 'precision','%1.f');

% dlmwrite('best_pairs_mst_corr_ers27_test2.txt', Best_pairs, 'delimiter', '\t', 'precision','%1.f');
return