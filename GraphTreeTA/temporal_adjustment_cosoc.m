%% temporal adjustment for Coso
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

printcleanfig = 1; % 1 = for papers, 0 = for research 
ylab = 'Volume Change'; yunits = 'm^3 \times 10^6'; scalefactor = 1.e6;

% data file name
%  datfile='mogi_vol_est.txt'
%datfile='cube_vol_est.txt'
% datfile='data_files/brick_vol_est_firstdraft.txt'
 datfile='brick_vol_est.txt'
 datfile = 'insar_time_series_data.csv';

dataset = datfile;

titlestring = sprintf('Estimates of Volume Change %s',datfile);




%% Read Data file
%  [VALS,NAMES] = read_table(datfile,1,4,0);
%  [ndat,ndummy] = size(VALS);  % number of pairs
% tm = cal2datetime(VALS(:,1));  % master epoch in decimal years
% ts = cal2datetime(VALS(:,2));  % slave epoch in decimal years
insar_table = readtable(datfile);
tm = cal2datetime(insar_table.caldate_master);  % master epoch in datetime
ts = cal2datetime(insar_table.caldate_slave);  % slave epoch in datetime
sat_id = zeros(numel(insar_table.bperp_m), 1);
sat_id(1:81) = 1;
sat_id(82:91) = 2;

% % truncate based on seis and prod data
% terupt = datetime(2005,01,01); % earliest seis/prod data
% iokm = find(tm > terupt);
% ioks = find(ts > terupt);
% iok = intersect(ioks,iokm);
% tm = tm(iok);
% ts = ts(iok);
% Qrate = Qrate(iok);
% rs = rs(iok);

ndat = numel(tm);

tm_cal = insar_table.caldate_master;
ts_cal = insar_table.caldate_slave;

tm_dt = tm;
tm_dt.Format = 'yyy-MM-dd';
ts_dt = ts;
ts_dt.Format = 'yyy-MM-dd';
tm = dyear(year(tm), month(tm), day(tm));
ts = dyear(year(ts), month(ts), day(ts));
% Qrate = VALS(:,3)./(ts-tm);
% rs = abs(VALS(:,4))./(ts-tm);
Qrate = insar_table.dVrate_m3peryr;
rs = insar_table.dVunc_m3peryr;


% %% test for dt vs dp
% 
% % set priors
% nu = 2.5E-01; % from best fit Mogi for stack
% G = 10e9; % shear modulus from best fit sun disk for stack
% E = 25e9; %G*2*(1+nu); 
% K = 5e9;%E/(3*(1-2*nu)); % assume biot coefficient of unity  => H = K (bulk modulus)
% mu_invH = 5e-10 %5e-10;% 5e-10;%1/K; %invPa
% sig_invH = .5*1/K; % invPa
% mu_pdot = -1.7237e+05; %-2.2e+05; %-9.3890e+05; %1.7237e+05;%-3.3e5; %Pa
% sig_pdot = (mu_pdot)/2; %Pa
% 
% mu_alphat = 1e-5; %invKel
% sig_alphat = .5*mu_alphat; % invKel
% mu_tdot = -1.4; %-2.6015; %degC
% sig_tdot = (mu_tdot)/2; %degC
% 
% % set volume
% %dx = 3.788e3 ; % best estimate by fitting cube to GIPhT best fit for mogi stack
% % dx = 3.1001780635E+03; % best fit cube for stack
% %V0 = dx^3;
% V0 = 3.0745993809E+03*2.2819266445E+03*1.9474147038E+03;%3.0745993809E+03*2.2819266445E+03*1.9474147038E+03; %6.8302e+09;
% 
% % results from TA 1 rate
% % TA_mu = -1.063458468050149e+06/V0;
% % TA_sig = 2.416898399574146e+04/V0;
% TA_mu = -1.2e6/V0%-1.076426787703105e+06/V0; %ENVI only
% TA_sig = 0.4e6/V0;%0.243102114744308e+04/V0;
% 
% % ENVI
% % test dT
% %alpha_t = 10^(-6); % from NASA proposal and Tabrez's modeling
% alpha_t = mu_alphat;
% test_dt = Qrate(ts < 2011)/(V0)/alpha_t; 
% test_dt_sig = abs(max(rs)/(V0)/mu_alphat);
% 
% figure;
% subplot(3, 2, 1)
% [at_mu, at_sig] = prodxypdf( mu_tdot, sig_tdot, mu_alphat, sig_alphat )
% histfit_with_prior_cosoc(test_dt*mu_alphat,ceil(sqrt(numel(test_dt))),'normal', at_mu, at_sig^2, (test_dt_sig*mu_alphat)^2, '1/yr', TA_mu, TA_sig, 0)
% xlabel('volumetric strain rate [1/yr]')
% xlabel('thermal volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated thermal volumetric strain rate (ENVI)')
% %printeps_e('hist_epsT_ENVI.eps')
% % figure
% % histfit_with_prior_cosoc(test_dt,ceil(sqrt(numel(test_dt))),'normal', mu_tdot, sig_tdot^2, (test_dt_sig)^2, 'degC/yr')
% % xlabel('temperature change rate [degC/yr]')
% % ylabel('count')
% % title('Estimated temperature change rate')
% % printeps_e('hist_dT.eps')
% 
% % test dP
% test_dp = Qrate(ts < 2011)/(V0)/(mu_invH);
% test_dp_sig = abs(max(rs)/(V0)/mu_invH);
% % figure;
% subplot(3, 2, 3)
% [hp_mu, hp_sig] = prodxypdf( mu_pdot, sig_pdot, mu_invH, sig_invH )
% histfit_with_prior_cosoc(test_dp*mu_invH,ceil(sqrt(numel(test_dp))),'normal', hp_mu, hp_sig^2, (test_dp_sig*mu_invH)^2, '1/yr', TA_mu, TA_sig, 0)
% xlabel('volumetric strain rate [1/yr]')
% xlabel('poroelastic volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated poroelastic volumetric strain rate (ENVI)')
% % printeps_e('hist_epsP_ENVI.eps')
% % figure
% % histfit_with_prior_cosoc(test_dp,ceil(sqrt(numel(test_dp))),'normal', mu_pdot, sig_pdot^2, (test_dp_sig)^2, 'Pa/yr', TA_mu/mu_invH, TA_sig/mu_invH)
% % xlabel('pressure change rate [Pa/yr]')
% % ylabel('count')
% % title('Estimated pressure change rate (ENVI)')
% % printeps_e('hist_dP_ENVI.eps')
% 
% % test both
% comb_mu = at_mu + hp_mu;
% comb_sig = sqrt(at_sig^2 + hp_sig^2);
% test_comb = Qrate(ts < 2011)/V0;
% test_comb_unc = max(rs)/V0;
% % figure; 
% subplot(3, 2, 5)
% histfit_with_prior_cosoc(test_comb,ceil(sqrt(numel(test_comb))),'normal', comb_mu, comb_sig^2, (test_comb_unc(1))^2, '1/yr', TA_mu, TA_sig, 2)
% xlabel('volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated volumetric strain rate (ENVI)')
% % printeps_e('hist_epsdTdP_ENVI.eps')
% 
% C = mu_alphat/mu_invH;
% dTdP = mu_tdot/mu_pdot;
% dTdPfac = C*dTdP;
% 
% dT_est = Qrate/V0/(mu_alphat+mu_invH*C/dTdPfac);%(mest)/(1+1/(dTdPfac*C));
% dP_est = dT_est*C/(dTdPfac);
% 
% eps_dT = mu_alphat*dT_est;
% eps_dP = mu_invH*dP_est;
% %dT_est = mest/(1+(mu_invH*mu_pdot)/(mu_alphat*mu_tdot));
% %dP_est = dT_est*(mu_invH*mu_pdot)/(mu_alphat*mu_tdot);
% 
% % figure
% % plot(eps_dP, eps_dT, 'k*')
% % hold on;
% % axislim = axis;
% % % lower bound for dP
% % line([hp_mu - hp_sig, hp_mu - hp_sig], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % upper bound for dP
% % line([hp_mu + hp_sig, hp_mu + hp_sig], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % lower bound for dT
% % line([axislim(1), axislim(2)], [at_mu - at_sig, at_mu - at_sig], 'LineStyle', '--')
% % % upper bound for dT
% % line([axislim(1), axislim(2)], [at_mu + at_sig, at_mu + at_sig], 'LineStyle', '--')
% % %ylim([1e-5, 8e-5])
% % plot(eps_dP, eps_dT, 'k*')
% % axis tight
% % xlabel('est. poroelastic vol. strain rate')
% % ylabel('est. thermal vol. strain rate')
% % printpdf('espdPvespdT_ENVI.pdf')
% % 
% % 
% % figure;
% % plot(dP_est, dT_est, 'k*')
% % hold on
% % axislim = axis;
% % % lower bound for dP
% % line([mu_pdot - sig_pdot, mu_pdot - sig_pdot], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % upper bound for dP
% % line([mu_pdot + sig_pdot, mu_pdot + sig_pdot], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % lower bound for dT
% % line([axislim(1), axislim(2)], [mu_tdot - sig_tdot, mu_tdot - sig_tdot], 'LineStyle', '--')
% % % upper bound for dT
% % line([axislim(1), axislim(2)], [mu_tdot + sig_tdot, mu_tdot + sig_tdot], 'LineStyle', '--')
% % plot(dP_est, dT_est, 'k*')
% % axis tight
% % xlabel('$\dot{P}_{est}$ [Pa/yr]', 'Interpreter', 'latex')
% % ylabel('$\dot{T}_{est}$ [degC/yr]', 'Interpreter', 'latex')
% % printpdf('dPvdT_ENVI.pdf')
% 
% % S1A
% % test dT
% %alpha_t = 10^(-6); % from NASA proposal and Tabrez's modeling
% V0 = 3.1119472417E+03*2.2886757881E+03*1.9536739341E+03;%3.0745993809E+03*2.2819266445E+03*1.9474147038E+03; %6.8302e+09;
% 
% 
% alpha_t = mu_alphat;
% test_dt = Qrate(tm >= 2011)/(V0)/alpha_t; 
% test_dt_sig = abs(max(rs)/(V0)/mu_alphat);
% 
% % results from TA 1 rate
% 
% TA_mu = -0.818329724891503e+06/V0; %S1A only
% TA_sig = 5e4/V0;%6.844960362538035e+04/V0;
% 
% 
% %%
% % figure;
% subplot(3, 2, 2)
% [at_mu, at_sig] = prodxypdf( mu_tdot, sig_tdot, mu_alphat, sig_alphat )
% histfit_with_prior_cosoc(test_dt*mu_alphat,3,'normal', at_mu, at_sig^2, (test_dt_sig*mu_alphat)^2, '1/yr', TA_mu, TA_sig, 0)
% 
% xlabel('volumetric strain rate [1/yr]')
% xlabel('thermal volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated thermal volumetric strain rate (S1A)')
% % printeps_e('hist_epsT_S1A.eps')
% % figure
% % histfit_with_prior_cosoc(test_dt,ceil(sqrt(numel(test_dt))),'normal', mu_tdot, sig_tdot^2, (test_dt_sig)^2, 'degC/yr')
% % xlabel('temperature change rate [degC/yr]')
% % ylabel('count')
% % title('Estimated temperature change rate')
% % printeps_e('hist_dT.eps')
% 
% % test dP
% test_dp = Qrate(tm >= 2011)/(V0)/(mu_invH);
% test_dp_sig = abs(max(rs)/(V0)/mu_invH);
% % figure;
% subplot(3, 2, 4)
% [hp_mu, hp_sig] = prodxypdf( mu_pdot, sig_pdot, mu_invH, sig_invH )
% histfit_with_prior_cosoc(test_dp*mu_invH,ceil(sqrt(numel(test_dp))),'normal', hp_mu, hp_sig^2, (test_dp_sig*mu_invH)^2, '1/yr', TA_mu, TA_sig, 0)
% xlabel('volumetric strain rate [1/yr]')
% xlabel('poroelastic volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated poroelastic volumetric strain rate (S1A)')
% printeps_e('hist_epsP_S1A.eps')
% % figure
% % histfit_with_prior_cosoc(test_dp,ceil(sqrt(numel(test_dp))),'normal', mu_pdot, sig_pdot^2, (test_dp_sig)^2, 'Pa/yr', TA_mu/mu_invH, TA_sig/mu_invH)
% % xlabel('pressure change rate [Pa/yr]')
% % ylabel('count')
% % title('Estimated pressure change rate')
% % printeps_e('hist_dP_S1A.eps')
% 
% % test both
% comb_mu = at_mu + hp_mu;
% comb_sig = sqrt(at_sig^2 + hp_sig^2);
% test_comb = Qrate(tm >= 2011)/V0;
% test_comb_unc = max(rs)/V0;
% % figure; 
% subplot(3, 2, 6)
% histfit_with_prior_cosoc(test_comb,ceil(sqrt(numel(test_comb))),'normal', comb_mu, comb_sig^2, (test_comb_unc(1))^2, '1/yr', TA_mu, TA_sig, 2)
% xlabel('volumetric strain rate [1/yr]')
% ylabel('count')
% title('Estimated volumetric strain rate (S1A)')
% % printeps_e('hist_epsdTdP_S1A.eps')
% %printeps_e('hist_eps_all.eps')
% 
% % C = mu_alphat/mu_invH;
% % dTdP = mu_tdot/mu_pdot;
% % dTdPfac = C*dTdP;
% % 
% % dT_est = Qrate(tm >= 2011)/V0/(mu_alphat+mu_invH*C/dTdPfac);%(mest)/(1+1/(dTdPfac*C));
% % dP_est = dT_est*C/(dTdPfac);
% % 
% % eps_dT = mu_alphat*dT_est;
% % eps_dP = mu_invH*dP_est;
% %dT_est = mest/(1+(mu_invH*mu_pdot)/(mu_alphat*mu_tdot));
% %dP_est = dT_est*(mu_invH*mu_pdot)/(mu_alphat*mu_tdot);
% % 
% % figure
% % plot(eps_dP, eps_dT, 'k*')
% % hold on;
% % axislim = axis;
% % % lower bound for dP
% % line([hp_mu - hp_sig, hp_mu - hp_sig], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % upper bound for dP
% % line([hp_mu + hp_sig, hp_mu + hp_sig], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % lower bound for dT
% % line([axislim(1), axislim(2)], [at_mu - at_sig, at_mu - at_sig], 'LineStyle', '--')
% % % upper bound for dT
% % line([axislim(1), axislim(2)], [at_mu + at_sig, at_mu + at_sig], 'LineStyle', '--')
% % %ylim([1e-5, 8e-5])
% % plot(eps_dP, eps_dT, 'k*')
% % axis tight
% % xlabel('est. poroelastic vol. strain rate')
% % ylabel('est. thermal vol. strain rate')
% % printpdf('espdPvespdT_S1A.pdf')
% % 
% % 
% % figure;
% % plot(dP_est, dT_est, 'k*')
% % hold on
% % axislim = axis;
% % % lower bound for dP
% % line([mu_pdot - sig_pdot, mu_pdot - sig_pdot], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % upper bound for dP
% % line([mu_pdot + sig_pdot, mu_pdot + sig_pdot], [axislim(3), axislim(4)], 'LineStyle', '--')
% % % lower bound for dT
% % line([axislim(1), axislim(2)], [mu_tdot - sig_tdot, mu_tdot - sig_tdot], 'LineStyle', '--')
% % % upper bound for dT
% % line([axislim(1), axislim(2)], [mu_tdot + sig_tdot, mu_tdot + sig_tdot], 'LineStyle', '--')
% % plot(dP_est, dT_est, 'k*')
% % axis tight
% % xlabel('$\dot{P}_{est}$ [Pa/yr]', 'Interpreter', 'latex')
% % ylabel('$\dot{T}_{est}$ [degC/yr]', 'Interpreter', 'latex')
% % printpdf('dPvdT_S1A.pdf')

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
% axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
hold on;
% ind23m = find(tm >= dyear(2016, 03, 14) & tm <=  dyear(2016, 03, 24));
% ind23s = find(ts >= dyear(2016, 03, 14) & ts <=  dyear(2016, 03, 24));
indother = 1:numel(tm);
% indother([ind23m; ind23s]) = [];

pest1 = mean(Qrate/scalefactor);
psig1 = std(Qrate/scalefactor);

% plot(tmid(indother(1)),Qrate(indother(1))/scalefactor, 'gs', 'MarkerFaceColor', 'g')
% plot(tmid(ind23s(1)),Qrate(ind23s(1))/scalefactor, 'r*')
% plot(tmid(ind23m(1)),Qrate(ind23m(1))/scalefactor, 'b+')
% plot((max(tu)+min(tu))/2.,pest1, 'ko', 'MarkerFaceColor', 'k')
% legend('normal operations pairs', 'pairs with slave in Stages2&3', 'pairs with master in Stages 2&3', 'mean value', 'Location', 'SouthWest')
errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
% errorbar_plus2(tmid(ind23s),Qrate(ind23s)/scalefactor,th(ind23s),rs(ind23s)/scalefactor,'r*',5);
% errorbar_plus2(tmid(ind23m),Qrate(ind23m)/scalefactor,th(ind23m),rs(ind23m)/scalefactor,'b+',5);


errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
% errorbar_plus2(tmid(ind23s),Qrate(ind23s)/scalefactor,th(ind23s),rs(ind23s)/scalefactor,'r*',5);
% errorbar_plus2(tmid(ind23m),Qrate(ind23m)/scalefactor,th(ind23m),rs(ind23m)/scalefactor,'b+',5);
ymin = nanmin(Qrate-rs)/scalefactor;
ymax = nanmax(Qrate+rs)/scalefactor;
% draw mean in blue
% color by stage dyear(2016, 03, 14), dyear(2016, 03, 24)

% [pest1,psig1,mse1] = lscov(ones(ndat,1),Qrate/scalefactor,diag(1./((rs/scalefactor).^2)));
% [pest1,psig1,mse1] = lscov((ts-tm),(VALS(:,3))/scalefactor,(1./(rs/scalefactor).^2));


%pest1 = mean(Qrate/scalefactor);
%psig1 = mean(Qrate/scalefactor);
errorbar_plus2((max(tu)+min(tu))/2.,pest1,(max(tu)-min(tu))/2.,psig1,'ko-',5);



% draw tbreaks in green
%tbreaks = [tu(1), dyear(2016, 02, 01), dyear(2016, 03, 14), dyear(2016, 03, 24), dyear(2016, 08, 23), tu(end)]; 
tbreaks = [tu(1), dyear(2017, 04, 27), tu(end)]; 

for i=1:numel(tbreaks)
    plot([tbreaks(i) tbreaks(i)],[ymin ymax],'k--');
end
hold off;
xlabel('Date');
ylabel(sprintf('Rate of %s [%s/year]',ylab,yunits));
set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
% datetick('x', 26);
if printcleanfig == 0
    title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
    
    printpdf(sprintf('%s_%sRATES.pdf',mfilename,dataset));
elseif printcleanfig == 1
    set(gca, 'FontSize', 14)
    %title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
    title('Estimated Rates of Volume Change for TSX pairs from 2017')
    [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
    tick_labels = datetime(axy, axm, axd);
%     tick_labels.Format = 'yyy-mmm-dd';
         set(gca, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
    printeps_e('Coso_mogi_RATES.eps');
    
end

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
dispflag = 1;

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
% Qdiff2 = adjustbp(tepochs,DD,Qdiff, trees,iuniqorbs, uniqdates);
Qdiff2 = adjustbp(tepochs,DD,insar_table.bperp_m, trees,iuniqorbs, uniqdates);
cal_date = [month(tepochs) day(tepochs)];
% plot 

tepochs = dyear(year(tepochs), month(tepochs), day(tepochs));

%h=plot_trees(tu, Qdiff2/scalefactor, DD, trees, 'calendar date [years]', sprintf('%s [%s]',ylab,yunits));

if printcleanfig == 0 %ktours = plotbp...
    plotbp(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), cal_date);
    printpdf(sprintf('%s_%s_TREES.pdf',mfilename,dataset),h);
elseif printcleanfig == 1
    figure;
    set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
    set(gcf,'DefaultTextInterpreter','tex');
   % plotbp_cosocMSF(tepochs, Qdiff2/scalefactor, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), cal_date);
    plotbp_cosocMSF(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), sat_id);
    printeps_e('Coso_mogi_TREES.eps');
end

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
tfuncs = {'nsegs'};
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
% tfuncs = {'nsegs'};
%tfuncs = {'exprdecay_coso2'};
%  tfuncs = {'logdecay'};
% tfuncs = {'sawtooth'};
% tfuncs = {'coso_exp2rate'};
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
        case {'exprdecay', 'logdecay'}
            metaparams(1) = tu(1);%dyear(1987, 01, 01);%tu(1); % reference time epoch
            %metaparams(1) = dyear(2003,1,1); % end of intrusion?
            metaparams(2) = inf;%6.5;%6.5;   % characteristic time constant in years
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);        
        case 'exprdecay_coso'
            metaparams(1) = tu(1); % reference time epoch
            %metaparams(1) = dyear(2003,1,1); % end of intrusion?
            metaparams(2) = inf;%6.5;%6.5;   % characteristic time constant in years
            metaparams(3) = inf;%6.5;%6.5;   % characteristic time constant in years
            metaparams(4) = V0;
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);        
       case 'coso_exp2rate'
           metaparams(1) = tu(1); % reference time epoch
           metaparams(2) = inf;   % characteristic time constant in years
          %metaparams(3) = 2003.6301 ; % first epoch in ENV data set
           %metaparams(3) = dyear(2003,1,1); % gap that is not spanned by data
           metaparams(3) = dyear(2011, 01, 01); % last epoch in ERS data set
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
        case 'exprdecay_coso2'
           metaparams(1) = tu(1); % reference time epoch
           metaparams(2) = inf;   % characteristic time constant in years
           metaparams(3) = dyear(2010, 01, 01); % tswitch
           metaparams(4) = inf; % end of intrusion?
         %   metaparams(4) = 2002.8630;%2003.6301 ; % first epoch in ENV data set
           metaparams(5) = tu(end); 
           
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
            tbreaks = [tu(1), dyear([2005:2011], 6*ones(size([2005:2011])), 1*ones(size([2005:2011]))), dyear(2014, 06, 01), dyear(2015, 06, 01), tu(end)];
            tbreaks = [tu(1), dyear([2005:2010], 6*ones(size([2005:2010])), 1*ones(size([2005:2010]))), tu(end)];
%                     tbreaks = [tu(1), dyear(2010, 01, 01), tu(end-7), tu(end-6), tu(end)];
            tbreaks = [tu(1), dyear(2010, 01, 01),  tu(end-7), tu(end-6), tu(end)];
 
        %   tbreaks = [tu(1), dyear(2010, 01, 01), tu(end)];
            tbreaks = unique(sort(tbreaks));
        case {'sawtooth'};
            metaparams = nan;
            load('Pumping/sawtooth_tbreaks.mat')
            tbreaks = unique(sort(tbreaks));
         %   tbreaks = sort([tu(1), dyear([2005:2015], 6*ones(size([2005:2012])), 1*ones(size([2005:2015]))), dyear([2005:2015], 12*ones(size([2005:2015])), 1*ones(size([2005:2015]))), tu(end)]);
            tbreaks = sort([tu(1), dyear([2005:2009], 6*ones(size([2005:2009])), 1*ones(size([2005:2009]))), dyear([2005:2009], 12*ones(size([2005:2009])), 1*ones(size([2005:2009]))), tu(end)]); % tu(end-7), tu(end-6), tu(end)]);
       %    tbreaks = sort([tu(1), dyear([2006:2010], 6*ones(size([2006:2010])), 1*ones(size([2006:2010]))), tu(end-7), tu(end-6), tu(end)]);

            % tbreaks = [tu(1), dyear(2010, 01, 01), dyear(2012, 01, 01), tu(end)];
        otherwise
            metaparams = nan;
            % few breaks
            tbreaks = [];
            tbreaks(end+1) = min(tu);
            tbreaks(end+1) = max(tu);
    end
   
    %% perform temporal adjustment
    tbreaks = colvec(sort(unique(tbreaks)))
    x_int = [10];
%     x_int = [10; mu_invH*5e4]; % for example, here we use a guess of 4 for tau and best switch time guess of June 6, 2002
%     x_int = [10; 10]; % for example, here we use a guess of 4 for tau and best switch time guess of June 6, 2002

    lb = [0];
    ub = [35];
%     lb = [10; (mu_invH - sig_invH)*1e4]; % column vector of lower bounds
%    ub = [30; (mu_invH + sig_invH)*1e5]; % column vector of upper bounds
%    lb = [0; 0]; % column vector of lower bounds
%    ub = [20; 20]; % column vector of upper bounds
% addpath('~/GraphTreeTA_UW_elena/Nonlin_TestCase/')
% [pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, SSWR, Vx, var, res_n, best_meta] = temporal_adjustment_nonlin(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams, x_int, lb, ub);

%    [pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, SSWR, Vx, var, res_n] = temporal_adjustment(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);
[pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, pfitd, rd, V, G, data, SSWR] = temporal_adjustment_function_coso(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);

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
    if strcmp(tfuncs, 'nsegs') == 1 | strcmp(tfuncs, 'sawtooth') == 1
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
    
    %% plot simple results
  %   save('insar_mod_pump.mat', 'pest', 'tbreaks', 'pfit', 'tfit', 'sigl', 'sigu')
  
        xlab = 'year'; 
   figure
%      load('TA_results_gps.mat')
% h =        plotpairs_coso2(tm_gps,ts_gps,Qdiff_gps*10/scalefactor,Qdsig_gps*10/scalefactor,Qdmod_gps*10/scalefactor ...
%              ,tfit_gps,pfit_gps*10/scalefactor,sigl_gps*10/scalefactor, sigu_gps*10/scalefactor ...
%              ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
%              ,'Cumulative Volume Change', tbreaks, 'mid');
%          hold on
h =       plotpairs_coso2(tm,ts,Qdiff/scalefactor,Qdsig/scalefactor,Qdmod/scalefactor ...
            ,tfit,pfit/scalefactor,sigl/scalefactor, sigu/scalefactor ...
            ,sprintf('%s', xlab),sprintf('%s [%s]',ylab,yunits)...
            ,'Cumulative volume produced', tbreaks, 'mid');
     printeps_e(sprintf('Coso_cuboid_%d%s',numel(pest), tfunc),h);    
end
   
diary off

return

