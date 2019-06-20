%% temporal adjustment for Dixie
% ECR 20170920

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
addpath('../Pumping')

printcleanfig = 1; % 1 = for papers, 0 = for research 
 ylab = 'Volume Change'; yunits = 'm^3 \times 10^6'; scalefactor = 1.e6;
%ylab = 'Volume Change'; yunits = 'm/yr'; scalefactor = 1;


% data file name
% datfile='mogi_vol_est.txt'
%datfile='mogi_vol_est_gps.txt'
datfile = 'brick_vol_est_gps_S1A.txt';
% datfile = 'brick_vol_est_gps_s1a_new.txt';
%datfile = 'brick_drho_est_gps.txt';
datfile = 'gps_post2010_time_series_data.csv';

dataset = datfile;

titlestring = sprintf('Estimates of Volume Change %s',datfile);




%% Read Data file
%  [VALS,NAMES] = read_table(datfile,0,4,0);
insar_table = readtable(datfile);

%  [ndat,ndummy] = size(VALS);  % number of pairs
% %  VALS = VALS(6:end, :);
% tm = cal2datetime(VALS(:,1));  % master epoch in decimal years
% ts = cal2datetime(VALS(:,2));  % slave epoch in decimal years
tm = cal2datetime(insar_table.caldate_master);  % master epoch in datetime
ts = cal2datetime(insar_table.caldate_slave);  % slave epoch in datetime


% % truncate based on seis and prod data
% terupt = datetime(2005,01,01); % earliest seis/prod data
% iokm = find(tm > terupt);
% ioks = find(ts > terupt);
% iok = intersect(ioks,iokm);
% tm = tm(iok);
% ts = ts(iok);
% Qrate = Qrate(iok);
% rs = rs(iok);
% tm_cal = VALS(:,1);
%  ts_cal = VALS(:,2);
tm_cal = insar_table.caldate_master;
ts_cal = insar_table.caldate_slave;
ndat = numel(tm);

tm_dt = tm;
tm_dt.Format = 'yyy-MM-dd';
ts_dt = ts;
ts_dt.Format = 'yyy-MM-dd';
tm = dyear(year(tm), month(tm), day(tm));
ts = dyear(year(ts), month(ts), day(ts));
%tm = VALS(:,1);
% %ts = VALS(:,2);
% Qrate = (VALS(:,3));%./(ts-tm);
% %rs = abs(max([abs(VALS(:,3)-VALS(:,4)), abs(VALS(:,3)-VALS(:,5))], [], 2))./(ts-tm);
% rs = sqrt(abs(VALS(:,4)))*1e3;%./(ts-tm);
Qrate = insar_table.dVrate_m3peryr;
rs = insar_table.dVunc_m3peryr;

% % test for dt vs dp
% dx = 3.1001780635E+03; % best fit cube for stack
% alpha_t = 10^(-5); % from NASA proposal and Tabrez's modeling
% test_dt = Qrate/(dx^3)/alpha_t; 
% figure; plot(test_dt)
% nu = 2.6927190156E-01; % from best fit Mogi for stack
% G = 4.423e9; % shear modulus from best fit sun disk for stack
% E = G*2*(1+nu); 
% K = E/(3*(1-2*nu)); % assume biot coefficient of unity  => H = K (bulk modulus)
% test_dp = Qrate/(dx^3)/(1/K);
% figure; plot(test_dp)

% unique epochs in years
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


% %% Print rates
% fprintf(1,'RATES:\ni,tm(i),ts(i),Qrate(i),rs(i)\n');
% for i=1:ndat
%      fprintf(1,'%3d %10.4f %10.4f %12.4e %12.4e\n',i,tm(i),ts(i),Qrate(i),rs(i));
% end

%% Plot rates
h_fig = figure; 
set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
set(h_fig,'DefaultTextInterpreter','tex');
% axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
hold on;
% ind23m = find(tm >= dyear(2016, 03, 14) & tm <=  dyear(2016, 03, 24));
% ind23s = find(ts >= dyear(2016, 03, 14) & ts <=  dyear(2016, 03, 24));
indother = 1:numel(tm);
% indother([ind23m; ind23s]) = [];

pest1 = mean(Qrate/scalefactor);
psig1 = std(Qrate/scalefactor);

% % plot(tmid(indother(1)),Qrate(indother(1))/scalefactor, 'gs', 'MarkerFaceColor', 'g')
% % plot(tmid(ind23s(1)),Qrate(ind23s(1))/scalefactor, 'r*')
% % plot(tmid(ind23m(1)),Qrate(ind23m(1))/scalefactor, 'b+')
% % plot((max(tu)+min(tu))/2.,pest1, 'ko', 'MarkerFaceColor', 'k')
% % legend('normal operations pairs', 'pairs with slave in Stages2&3', 'pairs with master in Stages 2&3', 'mean value', 'Location', 'SouthWest')
% errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
% % errorbar_plus2(tmid(ind23s),Qrate(ind23s)/scalefactor,th(ind23s),rs(ind23s)/scalefactor,'r*',5);
% % errorbar_plus2(tmid(ind23m),Qrate(ind23m)/scalefactor,th(ind23m),rs(ind23m)/scalefactor,'b+',5);
% 
% 
% errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
% % errorbar_plus2(tmid(ind23s),Qrate(ind23s)/scalefactor,th(ind23s),rs(ind23s)/scalefactor,'r*',5);
% % errorbar_plus2(tmid(ind23m),Qrate(ind23m)/scalefactor,th(ind23m),rs(ind23m)/scalefactor,'b+',5);
% ymin = nanmin(Qrate-rs)/scalefactor;
% ymax = nanmax(Qrate+rs)/scalefactor;
% % draw mean in blue
% % color by stage dyear(2016, 03, 14), dyear(2016, 03, 24)
% 
% % [pest1,psig1,mse1] = lscov(ones(ndat,1),Qrate/scalefactor,diag(1./((rs/scalefactor).^2)));
% % [pest1,psig1,mse1] = lscov((ts-tm),(VALS(:,3))/scalefactor,(1./(rs/scalefactor).^2));
% 
% 
% %pest1 = mean(Qrate/scalefactor);
% %psig1 = mean(Qrate/scalefactor);
% errorbar_plus2((max(tu)+min(tu))/2.,pest1,(max(tu)-min(tu))/2.,psig1,'ko-',5);
% 
% 
% 
% % draw tbreaks in green
% %tbreaks = [tu(1), dyear(2016, 02, 01), dyear(2016, 03, 14), dyear(2016, 03, 24), dyear(2016, 08, 23), tu(end)]; 
% tbreaks = [tu(1), dyear(2017, 04, 27), tu(end)]; 
% 
% for i=1:numel(tbreaks)
%     plot([tbreaks(i) tbreaks(i)],[ymin ymax],'k--');
% end
% hold off;
% xlabel('Date');
% ylabel(sprintf('Rate of %s [%s/year]',ylab,yunits));
% set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
% % datetick('x', 26);
% if printcleanfig == 0
%     title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
%     
%     printpdf(sprintf('%s_%sRATES.pdf',mfilename,dataset));
% elseif printcleanfig == 1
%     set(gca, 'FontSize', 14)
%     %title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
%     title('Estimated Rates of Volume Change for TSX pairs from 2017')
%     [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
%     tick_labels = datetime(axy, axm, axd);
% %     tick_labels.Format = 'yyy-mmm-dd';
%          set(gca, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
%     printeps_e('Coso_mogi_RATES.eps');
%     
% end

%% Convert rates into differences over time intervals
% this converts the volume from rate (m^3/year) to differential (m^3)
Qdiff = Qrate.* (ts - tm); % Qdiff is differential quantity observed over time interval
Qdsig = rs.* (ts - tm); % Qdsig is uncertainty (1-sigma) for differential observation

%% plot observations and costs
% figure; hold on;
% plot(nobs,cost,'r*');
% xlabel('number of observations');
% ylabel('cost');
% title(titlestring);
% printpdf(sprintf('%s_%sNOBSvsCOSTS.pdf',mfilename,dataset));


%% find and plot trees
dispflag = 0;

[tm_y, tm_m, tm_d] = dyear2date(tm);
[ts_y, ts_m, ts_d] = dyear2date(ts);
tm_dstruct = datetime(tm_y, tm_m, tm_d);
ts_dstruct = datetime(ts_y, ts_m, ts_d);

if exist('TA_gps_s1a.mat') == 2
    load('TA_gps_s1a.mat')
    [ndummy,mepochs] = size(DD);
    fprintf(1,'Number of unique epochs %d\n',mepochs);
    [ntrees,ndummy] = size(trees);
    fprintf(1,'Number of distinct trees %d\n',ntrees);
    [tepochs_yy, tepochs_mn, tepochs_day] = dyear2yyyymmdd(tepochs);
    cal_date = [tepochs_mn, tepochs_day];
else
%     [trees, DD tepochs, iepochs, iuniqorbs, uniqdates] = findtrees2(tm_dstruct',ts_dstruct');
    [trees, DD tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(tm, ts, dispflag,[1:numel(tu)-1] ,tm_cal,[2:numel(tu)],ts_cal);
    [ndummy,mepochs] = size(DD);
    fprintf(1,'Number of unique epochs %d\n',mepochs);
    [ntrees,ndummy] = size(trees);
    fprintf(1,'Number of distinct trees %d\n',ntrees);
    [te_y, te_m, te_d] = dyear2date(tepochs);
    %tm_dstruct = datetime(tm_y, tm_m, tm_d);
    cal_date = [te_m, te_d];
    %tepochs = dyear(year(tepochs), month(tepochs), day(tepochs));
    save('TA_gps_s1a.mat', 'DD', 'tepochs', 'trees')
end

%% naive temporal adjustment
%Qdiff2 = adjustbp(tepochs,DD,Qdiff, trees,iuniqorbs, uniqdates);
%cal_date = [month(tepochs) day(tepochs)];
%tepochs = dyear(year(tepochs), month(tepochs), day(tepochs));

%h=plot_trees(tu, Qdiff2/scalefactor, DD, trees, 'calendar date [years]', sprintf('%s [%s]',ylab,yunits));

% if printcleanfig == 0 %ktours = plotbp...
%     plotbp(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), cal_date);
%     printpdf(sprintf('%s_%s_TREES.pdf',mfilename,dataset),h);
% elseif printcleanfig == 1
%     figure;
%     set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
%     set(gcf,'DefaultTextInterpreter','tex');
%     plotbp(tepochs, Qdiff2/scalefactor, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), cal_date);
%     printeps_e('Coso_mogi_TREES.eps');
% end

% dispflag = 1;
% [trees, DD tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(tm,ts,dispflag);
% Qdiff2 = adjustbp(tepochs,DD,Qdiff, trees,iuniqorbs, uniqdates);
% plotbp(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits));
% printpdf(sprintf('%s_%s_TREES.pdf',mfilename,dataset));

%% loop over temporal fitting functions
%tfuncs = {'pwl','step-pwl','step-sec'};
%tfuncs = {'step'};
%tfuncs = {'step-pwl'};
%tfuncs = {'step-sec'};
%tfuncs = {'poly2','poly02'};
%tfuncs = {'nsegs','poly2','poly3','exp'};
%tfuncs = {'gaussian'};
%tfuncs = {'nsegs'};
%tfuncs = {'rate','step','pwl','step-pwl','step-sec','poly3','poly03'};
%tfuncs = {'rate','nsegs','step-pwl','exprdecay','okmokexp3'};
%  tfuncs = {'exprdecay'};
%tfuncs = {'pwl','poly2','poly3','exprdecay'};
%tfuncs = {'okmokexp3a'};
%tfuncs = {'okmokexp3'};
%tfuncs = {'rate','nsegs','exprdecay','okmokexp3'};
% tfuncs = {'poly2'};
%tfuncs = {'pwl', 'nsegs'};
%tfuncs = {'rate', 'nsegs', 'ber', 'ber_tikh', 'exprdecay','okmokexp4'};
 tfuncs = {'rate'};
%tfuncs = {'sawtooth'};
%tfuncs = {'rate', 'nsegs', 'pwl', 'exprdecay', 'okmokexp4'}; %full analysis
% dimension arrays
mses = nan(numel(tfuncs),1);
mparams = nan(numel(tfuncs),1);

for k = 1:numel(tfuncs)
    tfunc = tfuncs{k};
    diary(sprintf('%s_%s.out',mfilename,tfunc));
    
    fprintf(1,'\n-----------\n');
    
    fprintf(1,'Time function is now %s\n',tfunc);
    
   
    switch tfunc
        case 'poly2'
            metaparams(1) = terupt;
            tbreaks = tu(1:end);
        case 'exprdecay'
            metaparams(1) = terupt; % reference time epoch
            %metaparams(1) = dyear(2003,1,1); % end of intrusion?
            metaparams(2) = inf;%6.5;%6.5;   % characteristic time constant in years
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);           
       case 'okmokexp3'
           metaparams(1) = terupt; % reference time epoch
           metaparams(2) = 6.0;   % characteristic time constant in years
          %metaparams(3) = 2003.6301 ; % first epoch in ENV data set
           %metaparams(3) = dyear(2003,1,1); % gap that is not spanned by data
           metaparams(3) = 2002.8630; % last epoch in ERS data set
%          metaparams(3) = 2002.7151 % next to last epoch in ERS data set
%          metaparams(3) = 2002.61920 
          %metaparams(3) = dyear(2002,1,1); % gap spanned by at least one pair
          %metaparams(3) = 2003.0; % gap spanned by data
           
           metaparams(4) = 2.0;   % characteristic time constant in years
           %metaparams(5) = 7.e6; % [m^3] add arbitrary constant
           tbreaks = [];
           tbreaks(end+1) = min(tu);
           tbreaks(end+1) = max(tu);
           
            case 'okmokexp3a'
           metaparams(1) = terupt; % reference time epoch
           metaparams(2) = 4.5;   % characteristic time constant in years
          %metaparams(3) = 2003.6301 ; % first epoch in ENV data set
           %metaparams(3) = dyear(2003,1,1); % gap that is not spanned by data
           metaparams(3) = dyear(2002,5,1);%2002.8630; % last epoch in ERS data set
%          metaparams(3) = 2002.7151 % next to last epoch in ERS data set
%          metaparams(3) = 2002.61920 
          %metaparams(3) = dyear(2002,1,1); % gap spanned by at least one pair
          %metaparams(3) = 2003.0; % gap spanned by data


           metaparams(4) = 2.0;   % characteristic time constant in years
  
           metaparams(4) = 1.0;   % characteristic time constant in years
           %metaparams(5) = 7.e6; % [m^3] add arbitrary constant
           tbreaks = [];
           tbreaks(end+1) = min(tu);
           tbreaks(end+1) = max(tu);
        case 'okmokexp4'
           metaparams(1) = terupt; % reference time epoch
           metaparams(2) = 4.5;   % characteristic time constant in years
           metaparams(3) = dyear(2002,5,1); % beginning of intrusion?
           metaparams(4) = tu(end); % end of intrusion?
         %   metaparams(4) = 2002.8630;%2003.6301 ; % first epoch in ENV data set
           metaparams(5) = 4.5; 
           
           tbreaks = [];
           tbreaks(end+1) = min(tu);
           tbreaks(end+1) = max(tu);
        case 'ber_tikh'
            metaparams(1) = 1; % order for Tikhonov regularization
            metaparams(2) = 0.0000001; % lower bound for beta range
            metaparams(3) = 0.000001; % upper bound for beta range
            metaparams(4) = 0.00000001; % delta for beta range
            tbreaks = tu(1:end);
            case 'okmokexp5'
           metaparams(1) = terupt; % reference time epoch

           metaparams(2) = 2.0;   % characteristic time constant in years
           metaparams(3) = dyear(2002,1,1); % beginning of intrusion?
           metaparams(4) = 4;

           metaparams(2) = 6.0;   % characteristic time constant in years
           metaparams(3) = 2002.8630; % beginning of intrusion?
           metaparams(4) = 6;

          % metaparams(4) = dyear(2003,1,1); % end of intrusion?
%            metaparams(4) = 2003.6301 ; % first epoch in ENV data set

           tbreaks = [];
           tbreaks(end+1) = min(tu);
           tbreaks(end+1) = max(tu);
           
           case 'okmokexp6'
           metaparams(1) = terupt; % reference time epoch
           metaparams(2) = 2.0;   % characteristic time constant in years
           metaparams(3) = dyear(2002,1,1); % beginning of intrusion?
           metaparams(4) = tu(43);%dyear(2003,1,1); % end of intrusion?
%            metaparams(4) = 2003.6301 ; % first epoch in ENV data set

           tbreaks = [];
           tbreaks(end+1) = min(tu);
           tbreaks(end+1) = max(tu);
           
       case {'pwl', 'ber'};
            metaparams = nan;
            % many breaks
            tbreaks = [tu(1:end)];           
        case {'nsegs','step-pwl','step-sec'};
            metaparams = nan;
%             tbreaks = [tu(1), dyear(2016, 03, 14), dyear(2016, 03, 19), dyear(2016, 03, 24), tu(end)]; % 
            % tbreaks = [tu(1), dyear(2016, 03, 14), dyear(2016, 03, 24), tu(end)]; % 
            tbreaks = [tu(1), dyear([2005:2014], 6*ones(size([2005:2014])), 1*ones(size([2005:2014]))), tu(end)];
            tbreaks = unique(sort(tbreaks));
     case {'sawtooth'};
            metaparams = nan;       
            load('Pumping/sawtooth_tbreaks.mat')
            tbreaks = sort([tu(1), dyear([2011:2014], 6*ones(size([2011:2014])), 1*ones(size([2011:2014]))), dyear([2010:2014], 12*ones(size([2010:2014])), 1*ones(size([2010:2014]))), tu(end)]);
     %      tbreaks = sort([tu(1), dyear([2011:2014], 6*ones(size([2011:2014])), 1*ones(size([2011:2014]))), tu(end)]);

            %               tbreaks = sort([tu(1), dyear([2005:2008], 6*ones(size([2005:2008])), 1*ones(size([2005:2008]))), dyear([2005:2007], 12*ones(size([2005:2007])), 1*ones(size([2005:2007]))),...
         %                dyear([2011:2014], 6*ones(size([2011:2014])), 1*ones(size([2011:2014]))), dyear([2009:2014], 12*ones(size([2009:2014])), 1*ones(size([2009:2014]))), tu(end)]);

           % tbreaks = sort([tu(1), dyear([2005:2010], 6*ones(size([2005:2010])), 1*ones(size([2005:2010]))), dyear([2005:2009], 12*ones(size([2005:2009])), 1*ones(size([2005:2009]))), tu(end-7), tu(end-6), tu(end)]);
           
        otherwise
            metaparams = nan;
            % few breaks
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);
    end
   
    %% perform temporal adjustment
    tbreaks = colvec(sort(unique(tbreaks)))
    x_int = [4]; % for example, here we use a guess of 4 for tau and best switch time guess of June 6, 2002

    
    lb = [0]; % column vector of lower bounds
   ub = [10]; % column vector of upper bounds


%[pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, SSWR, Vx, var, res_n, best_meta] = temporal_adjustment_nonlin(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams, x_int, lb, ub);

   % [pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, SSWR, Vx, var, res_n] = temporal_adjustment(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);
[pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, pfitd, rd, V, G, data, SSWR] = temporal_adjustment_function_coso_gps_s1a(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);

   %   if strcmp(tfunc, 'pwl') ==1
%       pest_pwl=pest;
%   end
return
    mparam = numel(pest);
    %%
    % store statistics
    mses(k) = mse;
    mparams(k) = mparam;
    %tbreaks = [tu(1:end)];
    %% Calculate stats
    if strcmp(tfunc, 'nsegs') == 1
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

    end

      
xlab = 'year'; 
    figure;
    h  = plotpairs_cosoc_gps_s1a(tm,ts,Qdiff/scalefactor,Qdsig/scalefactor,Qdmod/scalefactor ...
            ,tfit,pfit/scalefactor,sigl/scalefactor, sigu/scalefactor ...
            ,sprintf('%s', xlab),sprintf('%s [%s]',ylab, yunits)...
            ,'Cumulative Volume Change', tbreaks, 'mid');

     printeps_e(sprintf('Coso_cuboid_%d%s_s1a',numel(pest), tfunc),h);    


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

diary off

return

