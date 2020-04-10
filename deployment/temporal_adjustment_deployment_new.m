%% temporal adjustment for Brady deployment
% last updated 20191013 Elena Reinisch

%% initialize
clear all; close all;

child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))
verbose = 1;
nf=0;

current_dir = pwd;
pdir = sprintf('%s/%s', pwd, 'Plots'); %name of folder to save prints to, if want in current working directory, replace with pwd
addpath(genpath(current_dir)); %calls on subdirectories as well (i.e. 'functions', etc)

printcleanfig = 1; % 1 = for papers, 0 = for research
ylab = 'Volume Change'; yunits = 'm^3 \times 10^4'; scalefactor = 1.e4;

% data file name
datfile='brick_vol_est.txt'

dataset = datfile;

titlestring = sprintf('Estimates of Volume Change %s',datfile);




%% Read Data file

insar_table = readtable('TSX_deployment_dV_ROI_sub_ext.txt', 'ReadVariableNames', 0);

tm_insar = table2array(insar_table(:,1));
ts_insar = table2array(insar_table(:,2));
dV_insar = table2array(insar_table(:,3));
dV_unc_insar = table2array(insar_table(:,4));

gps_table = readtable('brady_vol_est_gps_ROI_new_full2.txt', 'ReadVariableNames', 0);


% load GPS dates as dyear
tm_gps = table2array(gps_table(:,1));
ts_gps = table2array(gps_table(:,2));
% convert to Y M D
[tm_y, tm_m, tm_d] = dyear2yyyymmdd(tm_gps);
[ts_y, ts_m, ts_d] = dyear2yyyymmdd(ts_gps);
% convert to caldate
tm_gps = tm_y*1e4 + tm_m*1e2 +tm_d;
ts_gps = ts_y*1e4 + ts_m*1e2 +ts_d;

dV_gps = table2array(gps_table(:,3));
dV_unc_gps = table2array(gps_table(:,4));

%% join data sets
tm = [tm_insar; tm_gps];
ts = [ts_insar; ts_gps];
Qrate = [dV_insar; dV_gps];
rs = [dV_unc_insar; dV_unc_gps];

[ndat,ndummy] = size(tm);  % number of pairs
tm = cal2datetime(tm);  % master epoch in decimal years
ts = cal2datetime(ts);  % slave epoch in decimal years

ndat = numel(tm);

tm_dt = tm;
tm_dt.Format = 'yyy-MM-dd';
ts_dt = ts;
ts_dt.Format = 'yyy-MM-dd';
tm = dyear(year(tm), month(tm), day(tm));
ts = dyear(year(ts), month(ts), day(ts));

%% unique epochs in years
tu=sort(unique([tm ts]));
tu_dt = sort(unique([tm_dt ts_dt]));
tu_dt.Format = 'yyy-MM-dd';

tspan = sprintf('%10.4f to %10.4f\n',min(tu),max(tu));

% midpoints of each interval
tmid = (tm+ts)/2;

% half interval
th = (ts-tm)/2.;

% initial epoch
t0 = min(tu);

%% Print rates
fprintf(1,'RATES:\ni,tm(i),ts(i),Qrate(i),rs(i)\n');
for i=1:ndat
    fprintf(1,'%3d %10.4f %10.4f %12.4e %12.4e\n',i,tm(i),ts(i),Qrate(i),rs(i));
end

%% Plot rates
h_fig = figure;
set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
set(h_fig,'DefaultTextInterpreter','tex');

hold on;
indother = 1:numel(tm);

pest1 = mean(Qrate/scalefactor);
psig1 = std(Qrate/scalefactor);

errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
ymin = nanmin(Qrate-rs)/scalefactor;
ymax = nanmax(Qrate+rs)/scalefactor;

errorbar_plus2((max(tu)+min(tu))/2.,pest1,(max(tu)-min(tu))/2.,psig1,'ko-',5);

% draw tbreaks in green
tbreaks = [tu(1), dyear(2016, 02, 01), dyear(2016, 03, 14), dyear(2016, 03, 24), tu(end)];

for i=1:numel(tbreaks)
    plot([tbreaks(i) tbreaks(i)],[ymin ymax],'k--');
end
hold off;
xlabel('Date');
ylabel(sprintf('Rate of %s [%s/year]',ylab,yunits));
set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)

if printcleanfig == 0
    title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
    
    printpdf(sprintf('%s_%sRATES.pdf',mfilename,dataset));
elseif printcleanfig == 1
    set(gca, 'FontSize', 14)
    title('Estimated Rates of Volume Change for TSX and GPS during deployment')
    [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
    tick_labels = datetime(axy, axm, axd);
    set(gca, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
    printeps_e('brady_deployment_RATES.eps');
    
end

%% Convert rates into differences over time intervals
% this converts the volume from rate (m^3/year) to differential (m^3)
Qdiff = Qrate.* (ts - tm); % Qdiff is differential quantity observed over time interval
Qdsig = rs.* (ts - tm); % Qdsig is uncertainty (1-sigma) for differential observation

%% find and plot trees
dispflag = 0;

[tm_y, tm_m, tm_d] = dyear2date(tm);
[ts_y, ts_m, ts_d] = dyear2date(ts);
tm_dstruct = datetime(tm_y, tm_m, tm_d);
ts_dstruct = datetime(ts_y, ts_m, ts_d);

[trees, DD tepochs, iepochs, iuniqorbs, uniqdates] = findtrees2(tm_dstruct',ts_dstruct');

[ndummy,mepochs] = size(DD);
fprintf(1,'Number of unique epochs %d\n',mepochs);
[ntrees,ndummy] = size(trees);
fprintf(1,'Number of distinct trees %d\n',ntrees);

%% naive temporal adjustment
Qdiff2 = adjustbp(tepochs,DD,Qdiff, trees,iuniqorbs, uniqdates);
cal_date = [month(tepochs) day(tepochs)];
tepochs = dyear(year(tepochs), month(tepochs), day(tepochs));

%% loop over temporal fitting functions
tfuncs = {'nsegs'};

% dimension arrays
mses = nan(numel(tfuncs),1);
mparams = nan(numel(tfuncs),1);

for k = 1:numel(tfuncs)
    tfunc = tfuncs{k};
    diary(sprintf('%s_%s.out',mfilename,tfunc));
    
    fprintf(1,'\n-----------\n');
    
    fprintf(1,'Time function is now %s\n',tfunc);
    
    
    switch tfunc
        case {'nsegs'};
            metaparams = nan;
            tbreaks = [tu(1), dyear(2016, 03, 14), dyear(2016, 03, 19), dyear(2016, 03, 24), ...
                dyear(2018, 02, 22), dyear(2018, 05, 01), tu(end)];
            
            tbreaks = unique(sort(tbreaks));
        otherwise
            metaparams = nan;
            % few breaks
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);
    end
    
    %% perform temporal adjustment
    tbreaks = colvec(sort(unique(tbreaks)))
    [pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, SSWR, Vx, var, res_n] = temporal_adjustment(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);
    
    return
    
    mparam = numel(pest);
    %%
    % store statistics
    mses(k) = mse;
    mparams(k) = mparam;
    %tbreaks = [tu(1:end)];
    %% Calculate stats
    stats_mat = [];
    breaks = {};
    decision = {};
    for i = 2:numel(tbreaks)-1
        ind_1 = find(ts >= tbreaks(i-1) & ts < tbreaks(i));
        ind_2 = find(ts >= tbreaks(i) & ts < tbreaks(i+1));
        disp('number of data in first subset')
        numel(ind_1)
        disp('number of data in second subset')
        numel(ind_2)
        [ T_calc, df, T_alpha, Sp] = Ttest_2tail( pest(i-1), pest(i), 0, psig(i-1), psig(i), numel(ind_1), numel(ind_2), 0.05);
        % save stats
        if isnan(T_calc) == 1 || isnan(T_alpha) == 1
            H = 'test not valid';
        elseif abs(T_calc) < T_alpha
            H = 'fail to reject';
        else
            H = 'reject';
        end
        stats_mat(end+1,:) = [numel(ind_1), numel(ind_2), df, T_calc, T_alpha];
        breaks(end+1, :) = {sprintf('%7d', dyear2caldate(tbreaks(i-1))), sprintf('%7d', dyear2caldate(tbreaks(i)))},
        decision(end+1,1) = {H};
    end
    stats_table = table(breaks(:,1), breaks(:,2), stats_mat(:,1), stats_mat(:,2), stats_mat(:,3), stats_mat(:,4), stats_mat(:,5), decision, 'VariableNames', {'Break_1', 'Break_2', 'n_group1', 'n_group2', 'df', 'T_calc', 'T_alpha', 'H_decision'})
    stats_table_condensed =  table(breaks(:,1), breaks(:,2), stats_mat(:,3), stats_mat(:,4), stats_mat(:,5), decision, 'VariableNames', {'Break_1', 'Break_2',  'df', 'T_calc', 'T_alpha', 'H_decision'})
    
    return
    %% plot histogram of rate results by stage
    edges = -10:2.5:10;
    h_count = zeros(6, numel(edges)-1);
    nrates = zeros(6, 1);
    for i = 1:numel(tbreaks)-1
        stage_ind = i;
        nrates(i) = numel(Qrate(unique([find(tm >= tbreaks(stage_ind) & tm <= tbreaks(stage_ind+1))])));
        h_count(i,:) = histcounts(Qrate(unique([find(tm >= tbreaks(stage_ind) & tm <= tbreaks(stage_ind+1))]))/1e6, edges)./nrates(i);
        
    end
    figure;
    bar(edges(1:end-1),h_count');
    axis([-10, 10, 0, 1]);
    legend('normal operations (deployment stage 1)', 'site shutdown (deployment stage 2)', 'increased infield injection (deployment stage 3)', 'normal operations', 'site shutdown (2018)', 'normal operations');
    xlabel('rate of volume change [million m^3/yr]');
    ylabel('count');
    
    %% Plot Results multiple ways
    
    xlab = 'Date';
    if printcleanfig == 0
        %
        h  = plotpairs(tm,ts,Qdiff/scalefactor,Qdsig/scalefactor,Qdmod/scalefactor ...
            ,tfit,pfit/scalefactor,sigl/scalefactor, sigu/scalefactor ...
            ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
            ,strcat(titlestring,sprintf(' %s\nndat = %d mepochs = %d mparam = %d ntrees = %d rd = %d sqrt(MSE) = %8.2G'...
            ,tfunc, ndat, mepochs, mparam, ntrees, rd, sqrt(mse))), tbreaks, 'mast');
        
        
        % make PDF file
        printpdf(sprintf('%s_%s_%s_ADJUST.pdf',mfilename,dataset,tfunc),h);
    elseif printcleanfig == 1
        xlab = 'year';
        h  = plotpairs(tm,ts,Qdiff/scalefactor,Qdsig/scalefactor,Qdmod/scalefactor ...
            ,tfit,pfit/scalefactor,sigl/scalefactor, sigu/scalefactor ...
            ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
            ,'Cumulative Volume Change', tbreaks, 'mid');
        %             set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
        [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
        tick_labels = char(datetime(axy, axm, axd))
        set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
        ylim = [-5,5]
        %         ylabel('Volume Change [x 10^4 m^3]')
        %         set(gca,'DefaultTextInterpreter','latex');
        % make PDF file
        
        
        printeps_e(sprintf('brady_deployment_%d%s_seispumpcomp',numel(pest), tfunc),h);
        
    end
    
    %% plot TA with 2 plots
    ind1 = find(ts < dyear(2016, 07, 01));
    ind2 = find(tm >= dyear(2016, 07, 01));
    
    xlab = 'year';
    h  = plotpairs_gpsandinsar(mse, tm(ind1),ts(ind1),Qdiff(ind1)/scalefactor,Qdsig(ind1)/scalefactor,Qdmod(ind1)/scalefactor ...
        ,tfit(tfit < dyear(2016, 07, 01)),pfit(tfit < dyear(2016, 07, 01))/scalefactor,sigl(tfit < dyear(2016, 07, 01))/scalefactor, sigu(tfit < dyear(2016, 07, 01))/scalefactor ...
        ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
        ,'Cumulative Volume Change', tbreaks(1:4), 'mid');
    %             set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
    
    [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
    tick_labels = char(datetime(axy, axm, axd))
    set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    ylim = [-5,5]
    %         ylabel('Volume Change [x 10^4 m^3]')
    %         set(gca,'DefaultTextInterpreter','latex');
    % make PDF file
    %
    
    tfit2 = find(tfit >= dyear(2016, 07, 01));
    h  = plotpairs_gpsandinsar(mse, tm(ind2),ts(ind2),Qdiff(ind2)/scalefactor,Qdsig(ind2)/scalefactor,Qdmod(ind2)/scalefactor - Qdmod(ind2(1))/scalefactor ...
        ,tfit(tfit >= dyear(2016, 07, 01)),pfit(tfit2)/scalefactor - pfit(tfit2(1))/scalefactor,(sigl(tfit2) - sigl(tfit2(1)))/scalefactor, (sigu(tfit2) - sigu(tfit2(1)))/scalefactor ...
        ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
        ,'Cumulative Volume Change from 1 July 2016', tbreaks(5:end), 'mid');
    %             set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
    axis 'tight'
    [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
    % return
    tick_labels = char(datetime(axy, axm, axd))
    set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    
    ylim = [-5,5]
    
    
    
    %% plot rates
    figure;
    
    tmid_insar = datetime2dyear(cal2datetime(tm_insar)) + (datetime2dyear(cal2datetime(ts_insar)) - datetime2dyear(cal2datetime(tm_insar)))/2;
    tmid_gps = datetime2dyear(cal2datetime(tm_gps)) + (datetime2dyear(cal2datetime(ts_gps)) - datetime2dyear(cal2datetime(tm_gps)))/2;
    %
    % Qrate = [dV_insar; dV_gps];%./(ts-tm);
    % % rs = 8*[dV_unc_insar; dV_unc_gps*1e-4];%*.25e-3];%./(ts-tm);
    % rs = [dV_unc_insar; dV_unc_gps];%*.25e-3];%./(ts-tm);
    
    sf = 1e5;
    
    figure;
    hold on
    set(gca, 'FontSize', 12)
    
    % plot gps and insar
    %line([tmid_insar, tmid_insar]', [dV_insar-dV_unc_insar, dV_insar+dV_unc_insar]'/sf, 'Color', 'k')
    %line([tmid_gps, tmid_gps]', [dV_gps-sqrt(dV_unc_gps), dV_gps+sqrt(dV_unc_gps)]'/sf, 'Color', 'b')
    pair1_ind = find(tmid_insar <= tbreaks(2) & tmid_insar >= tbreaks(1));
    pair2_ind = find(tmid_insar <= tbreaks(3) & tmid_insar >= tbreaks(2));
    pair3_ind = find(tmid_insar <= tbreaks(4) & tmid_insar >= tbreaks(3));
    pair4_ind = find(tmid_insar <= tbreaks(5) & tmid_insar >= tbreaks(4));
    pair5_ind = find(tmid_insar <= tbreaks(6) & tmid_insar >= tbreaks(5));
    pair6_ind = find(tmid_insar <= tbreaks(7) & tmid_insar >= tbreaks(6));
    pair1_ind = find(tmid_gps <= tbreaks(2) & tmid_gps >= tbreaks(1));
    pair2_ind = find(tmid_gps <= tbreaks(3) & tmid_gps >= tbreaks(2));
    pair3_ind = find(tmid_gps <= tbreaks(4) & tmid_gps >= tbreaks(3));
    pair4_ind = find(tmid_gps <= tbreaks(5) & tmid_gps >= tbreaks(4));
    pair5_ind = find(tmid_gps <= tbreaks(6) & tmid_gps >= tbreaks(5));
    pair6_ind = find(tmid_gps <= tbreaks(7) & tmid_gps >= tbreaks(6));
    
    %     for i = 1:numel(tbreaks)-1
    %        pair_ind = find(tmid_gps < tbreaks(i+1) & tmid_gps > tbreaks(i));
    %        boxplot(dV_gps(pair_ind)/sf, 'positions', tbreaks(i)+(tbreaks(i+1) - tbreaks(i))/2)
    %        axis([dyear(2016, 3, 1), dyear(2018, 10, 2), -200, 200])
    %     end
    %
    %    return
    errorbar_plus2(tbreaks(1)+(tbreaks(2) - tbreaks(1))/2, pest(1)/sf, (tbreaks(2) - tbreaks(1))/2, psig(1)/sf, 'go', 6, 3)
    errorbar_plus2(tbreaks(2)+(tbreaks(3) - tbreaks(2))/2, pest(2)/sf, (tbreaks(3) - tbreaks(2))/2, psig(2)/sf, 'ro', 6, 3)
    errorbar_plus2(tbreaks(3)+(tbreaks(4) - tbreaks(3))/2, pest(3)/sf, (tbreaks(4) - tbreaks(3))/2, psig(3)/sf, 'bo', 6, 3)
    errorbar_plus2(tbreaks(4)+(tbreaks(2) - tbreaks(1))/2, pest(4)/sf, (tbreaks(2) - tbreaks(1))/2, psig(4)/sf, 'go', 6, 3)
    errorbar_plus2(tbreaks(5)+(tbreaks(6) - tbreaks(5))/2, pest(5)/sf, (tbreaks(6) - tbreaks(5))/2, psig(2)/sf, 'ro', 6, 3)
    errorbar_plus2(tbreaks(6)+(tbreaks(7) - tbreaks(6))/2, pest(6)/sf, (tbreaks(7) - tbreaks(6))/2, psig(4)/sf, 'go', 6, 3)
    hold on
    for i=1:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[-15, 15],'k--');
    end
    axis tight
    
    set(gca, 'XTick', tbreaks)
    [axy, axm, axd] = dyear2yyyymmdd(str2double(get(gca, 'XTickLabel')));
    tick_labels = datetime(axy, axm, axd);
    tick_labels.Format = 'yy-MM-dd';
    set(gca, 'XTickLabel',datestr(tick_labels), 'XTickLabelRotation', -75) %char(tick_labels)
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    
    ylabel('Volume Change Rate [x 10^5 m^3/yr]')
    %return
    subtightplot(1, 2, 2)
    hold on
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    errorbar_plus2(tbreaks(5)-(tbreaks(5) - tbreaks(4))/6, pest(4)/sf, (tbreaks(5) - tbreaks(4))/6, psig(4)/sf, 'go', 6, 3)
    if numel(tbreaks) > 5
        errorbar_plus2(tbreaks(5)+(tbreaks(6) - tbreaks(5))/2, pest(5)/sf, (tbreaks(6) - tbreaks(5))/2, psig(2)/sf, 'ro', 6, 3)
        errorbar_plus2(tbreaks(6)+(tbreaks(7) - tbreaks(6))/2, pest(6)/sf, (tbreaks(7) - tbreaks(6))/2, psig(4)/sf, 'go', 6, 3)
    end
    %     ylim([-8,10]);
    for i=5:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[-15, 15],'k--');
    end
    %         for i=1:numel(tbreaks)-1
    %             plot([tbreaks(i) tbreaks(i)],[-8, 10],'k--');
    %         end
    %         plot([tbreaks(6)+(tbreaks(7) - tbreaks(6))/4 tbreaks(6)+(tbreaks(7) - tbreaks(6))/4],[-8, 10],'k--');
    axis tight
    %    set(gca, 'YTick', 'None')
    set(gca, 'YTickLabel', '{}')
    set(gca, 'XTick', tbreaks(5:end))
    [axy, axm, axd] = dyear2yyyymmdd(str2double(get(gca, 'XTickLabel')));
    tick_labels = datetime(axy, axm, axd);
    tick_labels.Format = 'yy-MM-dd';
    set(gca, 'XTickLabel',datestr(tick_labels), 'XTickLabelRotation', -75) %char(tick_labels)
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    % set(gca, 'XTick', [tbreaks(1:end-1), tbreaks(6)+(tbreaks(7) - tbreaks(6))/4])
    % return
    cax = breakxaxis([dyear(2016,04,16), dyear(2018,02,01)], 0.0015)
    
    
    %title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
    %title(sprintf('Estimated Rates of Volume Change \n for TSX and GPS during deployment'))
    ylabel('Volume Change Rate [x 10^5 m^3/yr]')
    
    %     [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTickLabel'));
    %     tick_labels = datetime(axy, axm, axd);
    [axy, axm, axd] = dyear2yyyymmdd(cax.leftAxes.XTick);
    tick_labels = datetime(axy, axm, axd);
    cax.leftAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    cax.leftAxes.XTickLabelRotation = -50;
    [axy, axm, axd] = dyear2yyyymmdd(cax.rightAxes.XTick);
    tick_labels = datetime(axy, axm, axd);
    
    cax.rightAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    cax.rightAxes.XTickLabelRotation = -50;
    %       cax2 = breakxaxis([dyear(2018,05,22), dyear(2018,09,11)], 0.0015)
    %        [axy, axm, axd] = dyear2yyyymmdd(cax2.leftAxes.XTick);
    %     tick_labels = datetime(axy, axm, axd);
    %     cax2.leftAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    %     cax2.leftAxes.XTickLabelRotation = -50;
    %     [axy, axm, axd] = dyear2yyyymmdd(cax2.rightAxes.XTick);
    %     tick_labels = datetime(axy, axm, axd);
    %     cax2.rightAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    %     cax2.rightAxes.XTickLabelRotation = -50;
    
    
    
    %     tick_labels.Format = 'yyy-mmm-dd';
    %     set(cax, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
    printeps_e('brady_deployment_RATES_2.eps');
    
    
    
    
    
    
    %     set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    %        set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    %     errorbar_plus2(tbreaks(1:end-1)+diff(tbreaks)/2, pest, diff(tbreaks)/2, psig, 'bo')
    sf = 1e5;
    % %     plot(tbreaks(1)+(tbreaks(2) - tbreaks(1))/2, pest(1)/sf,'g-', 'LineWidth', 2);
    % %     hold on
    % %     plot(tbreaks(2)+(tbreaks(3) - tbreaks(2))/2, pest(2)/sf,'r-', 'LineWidth', 2);
    % %     plot(tbreaks(3)+(tbreaks(4) - tbreaks(3))/2, pest(3)/sf,'b-', 'LineWidth', 2);
    % %      plot(tbreaks(4)+(tbreaks(5) - tbreaks(4))/2, pest(4)/sf,'g-', 'LineWidth', 2);
    %
    %     errorbar_plus2(tbreaks(1)+(tbreaks(2) - tbreaks(1))/2, pest(1)/sf, (tbreaks(2) - tbreaks(1))/2, psig(1)/sf, 'go', 6, 3)
    %     errorbar_plus2(tbreaks(2)+(tbreaks(3) - tbreaks(2))/2, pest(2)/sf, (tbreaks(3) - tbreaks(2))/2, psig(2)/sf, 'ro', 6, 3)
    %     errorbar_plus2(tbreaks(3)+(tbreaks(4) - tbreaks(3))/2, pest(3)/sf, (tbreaks(4) - tbreaks(3))/2, psig(3)/sf, 'bo', 6, 3)
    %     errorbar_plus2(tbreaks(4)+(tbreaks(5) - tbreaks(4))/2, pest(4)/sf, (tbreaks(5) - tbreaks(4))/2, psig(4)/sf, 'go', 6, 3)
    % %
    % %     legend('normal operations', 'site shutdown', sprintf('increased injection \n & pulsing'))
    %
    % return
    
    
    pair1_ind = find(tmid_insar <= tbreaks(2) & tmid_insar >= tbreaks(1));
    pair2_ind = find(tmid_insar <= tbreaks(3) & tmid_insar >= tbreaks(2));
    pair3_ind = find(tmid_insar <= tbreaks(4) & tmid_insar >= tbreaks(3));
    pair4_ind = find(tmid_insar <= tbreaks(5) & tmid_insar >= tbreaks(4));
    pair5_ind = find(tmid_insar <= tbreaks(6) & tmid_insar >= tbreaks(5));
    pair6_ind = find(tmid_insar <= tbreaks(7) & tmid_insar >= tbreaks(6));
    pair1_ind = find(tmid_gps <= tbreaks(2) & tmid_gps >= tbreaks(1));
    pair2_ind = find(tmid_gps <= tbreaks(3) & tmid_gps >= tbreaks(2));
    pair3_ind = find(tmid_gps <= tbreaks(4) & tmid_gps >= tbreaks(3));
    pair4_ind = find(tmid_gps <= tbreaks(5) & tmid_gps >= tbreaks(4));
    pair5_ind = find(tmid_gps <= tbreaks(6) & tmid_gps >= tbreaks(5));
    pair6_ind = find(tmid_gps <= tbreaks(7) & tmid_gps >= tbreaks(6));
    
    for i = 1:numel(tbreaks)-1
        pair_ind = find(tmid_gps < tbreaks(i+1) & tmid_gps > tbreaks(i));
        boxplot(dV_gps(pair_ind)/sf, 'positions', tbreaks(i)+(tbreaks(i+1) - tbreaks(i))/2)
        axis([dyear(2016, 3, 1), dyear(2018, 10, 2), -200, 200])
    end
    
    
    %%
    subtightplot(1, 2, 1)
    
    hold on
    set(gca, 'FontSize', 12)
    
    pair1_ind = find(tmid_insar <= tbreaks(2) & tmid_insar >= tbreaks(1));
    pair2_ind = find(tmid_insar <= tbreaks(3) & tmid_insar >= tbreaks(2));
    pair3_ind = find(tmid_insar <= tbreaks(4) & tmid_insar >= tbreaks(3));
    pair4_ind = find(tmid_insar <= tbreaks(5) & tmid_insar >= tbreaks(4));
    pair1g_ind = find(tmid_gps <= tbreaks(2) & tmid_gps >= tbreaks(1));
    pair2g_ind = find(tmid_gps <= tbreaks(3) & tmid_gps >= tbreaks(2));
    pair3g_ind = find(tmid_gps <= tbreaks(4) & tmid_gps >= tbreaks(3));
    pair4g_ind = find(tmid_gps <= tbreaks(5) & tmid_gps >= tbreaks(4));
    %     for i = 1:4-1
    %         axis([dyear(2016, 3, 1), tbreaks(4)+(tbreaks(2) - tbreaks(1)), -200, 200])
    %        pair_ind = find(tmid_gps < tbreaks(i+1) & tmid_gps > tbreaks(i));
    %        boxplot(dV_gps(pair_ind)/sf, 'position', i, 'boxstyle', 'outline', 'width', .05)
    %       % axis([dyear(2016, 3, 1), tbreaks(4)+(tbreaks(2) - tbreaks(1)), -200, 200])
    %     end
    
    pair_indg = [pair1g_ind; pair2g_ind; pair3g_ind; pair4g_ind];
    gg = [(tbreaks(1)+(tbreaks(2) - tbreaks(1))/2)*ones(length(pair1g_ind), 1); (tbreaks(2)+(tbreaks(3) - tbreaks(2))/2)*ones(numel(pair2g_ind), 1);...
        (tbreaks(3)+(tbreaks(4) - tbreaks(3))/2)*ones(numel(pair3g_ind), 1); (tbreaks(4)+(tbreaks(2) - tbreaks(1))/2)*ones(numel(pair4g_ind), 1)];
    boxplot(dV_gps(pair_indg)/sf, gg, 'plotstyle', 'compact', 'color', [0.5, 0.5, 0.5], 'positions', unique(gg), 'boxstyle', 'outline', 'widths', 0.005, 'symbol', '+')
    hAx = gca;
    hAx.YLim = [-90, 90];
    hAx.XTick = unique(gg);
    xtk=hAx.XTick;
    hold on
    pair_ind = [pair1_ind; pair2_ind; pair3_ind; pair4_ind];
    gi = [(tbreaks(1)+(tbreaks(2) - tbreaks(1))/2)*ones(length(pair1_ind), 1); (tbreaks(2)+(tbreaks(3) - tbreaks(2))/2)*ones(numel(pair2_ind), 1);...
        (tbreaks(3)+(tbreaks(4) - tbreaks(3))/2)*ones(numel(pair3_ind), 1); (tbreaks(4)+(tbreaks(2) - tbreaks(1))/2)*ones(numel(pair4_ind), 1)];
    hi = boxplot(dV_insar(pair_ind)/sf, gi, 'plotstyle', 'compact', 'color', [0.4940, 0.1840, 0.5560], 'positions', unique(gi),'boxstyle', 'outline', 'widths', 0.005, 'symbol', 'd')
    set(hi,'LineWidth',2)
    hAx = gca;
    hAx.YLim = [-90, 90];
    hAx.XTick = unique(gg);
    xtk=hAx.XTick;
    hAx.XLim = [xtk(1) - (tbreaks(2) - tbreaks(1))/2, xtk(4) + (tbreaks(2) - tbreaks(1))/2];
    hold on
    errorbar_plus2(xtk(1), pest(1)/sf, (tbreaks(2) - tbreaks(1))/2, psig(1)/sf, 'g', 6, 3)
    errorbar_plus2(xtk(2), pest(2)/sf, (tbreaks(3) - tbreaks(2))/2, psig(2)/sf, 'r', 6, 3)
    errorbar_plus2(xtk(3), pest(3)/sf, (tbreaks(4) - tbreaks(3))/2, psig(3)/sf, 'b', 6, 3)
    errorbar_plus2(xtk(4), pest(4)/sf, (tbreaks(2) - tbreaks(1))/2, psig(4)/sf, 'g', 6, 3)
    
    
    
    
    
    
    %     ylim([-8,10]);
    hold on
    for i=1:4
        plot([tbreaks(i) tbreaks(i)],[-90, 90],'k:', 'LineWidth', .1);
    end
    
    axis tight
    
    set(gca, 'XTick', tbreaks(1:4))
    
    %      tick_labels = datetime([2016, 2016, 2016, 2016], [3, 3, 3, 3], [1, 14, 19, 24]);
    %      tick_labels.Format = 'yyyy-mm-dd';
    tick_labels = {'2016-Mar-01', '2016-Mar-14', '2016-Mar-19', '2016-Mar-24'};
    set(gca, 'XTickLabel',tick_labels, 'XTickLabelRotation', -75) %char(tick_labels)
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('Volume Change Rate [x 10^5 m^3/yr]')
    
    
    
    subtightplot(1, 2, 2)
    hold on
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    
    
    pair5_ind = find(tmid_insar <= tbreaks(6) & tmid_insar >= tbreaks(5));
    pair6_ind = find(tmid_insar <= tbreaks(7) & tmid_insar >= tbreaks(6));
    pair5g_ind = find(tmid_gps <= tbreaks(6) & tmid_gps >= tbreaks(5));
    pair6g_ind = find(tmid_gps <= tbreaks(7) & tmid_gps >= tbreaks(6));
    
    pair_ind = [pair5_ind; pair6_ind];
    gi = [(tbreaks(5)+(tbreaks(6) - tbreaks(5))/2)*ones(numel(pair5_ind), 1);...
        (tbreaks(6)+(tbreaks(7) - tbreaks(6))/2)*ones(numel(pair6_ind), 1)];
    hi = boxplot(dV_insar(pair_ind)/sf, gi, 'plotstyle', 'compact', 'color', [0.4940, 0.1840, 0.5560], 'positions', unique(gi),'boxstyle', 'outline', 'widths', 0.05, 'symbol', 'd')
    set(hi,'LineWidth',2)
    hAx = gca;
    hAx.YLim = [-90, 90];
    hAx.XTick = unique(gi);
    xtk=hAx.XTick;
    hold on
    
    pair_indg = [pair5g_ind; pair6g_ind];
    gg = [(tbreaks(5)+(tbreaks(6) - tbreaks(5))/2)*ones(numel(pair5g_ind), 1);...
        (tbreaks(6)+(tbreaks(7) - tbreaks(6))/2)*ones(numel(pair6g_ind), 1)];
    boxplot(dV_gps(pair_indg)/sf, gg, 'plotstyle', 'compact', 'color', [0.5, 0.5, 0.5], 'positions', unique(gg), 'boxstyle', 'outline', 'widths', 0.05, 'symbol', '+')
    %    boxplot(dV_insar(pair_ind)/sf, gi, 'plotstyle', 'compact', 'color', [0.4940, 0.1840, 0.5560], 'positions', unique(gi),'boxstyle', 'outline', 'widths', 0.05, 'symbol', 'd')
    
    
    hAx = gca;
    hAx.YLim = [-90, 90];
    hAx.XTick = unique(gi);
    xtk=hAx.XTick;
    hAx.XLim = [xtk(1) - (tbreaks(5) - tbreaks(4))/6, xtk(2) + (tbreaks(7) - tbreaks(6))/2];
    hold on
    
    
    
    errorbar_plus2(tbreaks(5)-(tbreaks(5) - tbreaks(4))/80, pest(4)/sf, (tbreaks(5) - tbreaks(4))/80, psig(4)/sf, 'g', 6, 3)
    if numel(tbreaks) > 5
        errorbar_plus2(tbreaks(5)+(tbreaks(6) - tbreaks(5))/2, pest(5)/sf, (tbreaks(6) - tbreaks(5))/2, psig(2)/sf, 'r', 6, 3)
        errorbar_plus2(tbreaks(6)+(tbreaks(7) - tbreaks(6))/2, pest(6)/sf, (tbreaks(7) - tbreaks(6))/2, psig(4)/sf, 'g', 6, 3)
    end
    %     ylim([-8,10]);
    for i=5:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[-90, 90],'k:');
    end
    
    axis tight
    set(gca, 'YTickLabel', '{}')
    set(gca, 'YTick', [])
    %    set(gca, 'XTick', tbreaks(5:end))
    %    [axy, axm, axd] = dyear2yyyymmdd(str2double(get(gca, 'XTickLabel')));
    %     tick_labels = datetime(axy, axm, axd);
    %     tick_labels.Format = 'yy-MM-dd';
    % 	set(gca, 'XTickLabel',datestr(tick_labels), 'XTickLabelRotation', -75) %char(tick_labels)
    %  set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    %        set(gca, 'XTickLabel', datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75)
    %    % set(gca, 'XTick', [tbreaks(1:end-1), tbreaks(6)+(tbreaks(7) - tbreaks(6))/4])
    set(gca, 'XTick', tbreaks(5:7))
    
    %      tick_labels = datetime([2016, 2016, 2016, 2016], [3, 3, 3, 3], [1, 14, 19, 24]);
    %      tick_labels.Format = 'yyyy-mm-dd';
    tick_labels = {'2018-Feb-21', '2018-Apr-30', '2018-Oct-02'};
    set(gca, 'XTickLabel',tick_labels, 'XTickLabelRotation', -75) %char(tick_labels)
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    %     ylabel('Volume Change Rate [x 10^5 m^3/yr]')
    
    
    return
    cax = breakxaxis([dyear(2016,04,16), dyear(2018,02,01)], 0.0015)
    
    
    %title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
    %title(sprintf('Estimated Rates of Volume Change \n for TSX and GPS during deployment'))
    ylabel('Volume Change Rate [x 10^5 m^3/yr]')
    
    %     [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTickLabel'));
    %     tick_labels = datetime(axy, axm, axd);
    [axy, axm, axd] = dyear2yyyymmdd(cax.leftAxes.XTick);
    tick_labels = datetime(axy, axm, axd);
    cax.leftAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    cax.leftAxes.XTickLabelRotation = -50;
    [axy, axm, axd] = dyear2yyyymmdd(cax.rightAxes.XTick);
    tick_labels = datetime(axy, axm, axd);
    
    cax.rightAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    cax.rightAxes.XTickLabelRotation = -50;
    %       cax2 = breakxaxis([dyear(2018,05,22), dyear(2018,09,11)], 0.0015)
    %        [axy, axm, axd] = dyear2yyyymmdd(cax2.leftAxes.XTick);
    %     tick_labels = datetime(axy, axm, axd);
    %     cax2.leftAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    %     cax2.leftAxes.XTickLabelRotation = -50;
    %     [axy, axm, axd] = dyear2yyyymmdd(cax2.rightAxes.XTick);
    %     tick_labels = datetime(axy, axm, axd);
    %     cax2.rightAxes.XTickLabel = datestr(tick_labels, 'yy-mmm-dd');
    %     cax2.rightAxes.XTickLabelRotation = -50;
    
    
    
    % %     tick_labels.Format = 'yyy-mmm-dd';
    % %     set(cax, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
    %     printeps_e('brady_deployment_RATES_2.eps');
    
    
    %% display results
    [nd, md] = size(DD);
    fprintf(1,'PARAMETERS:\ni,pest(i),psig(i)\n');
    for i=1:mparam
        fprintf(1,'%3d %10.4f %12.4e\n',i,pest(i),psig(i));
    end
    
end

% print summary
fprintf(1,'\n\n------\nSUMMARY:\n');
fprintf(1,'Number of data %d\n',ndat);
fprintf(1,'Number of unique epochs %d\n',mepochs);
fprintf(1,'Number of trees (distinct trees) %d\n',ntrees);
fprintf(1,'K, time function, mparam, sqrt(MSE)\n');
for k=1:numel(tfuncs)
    fprintf(1,'%d %#20s %3d %10.4f\n',k,tfuncs{k},mparams(k),sqrt(mses(k)));
end
res = (Qdiff - G*pest);
% save('test.mat', 'res')
diary off

return

