function h_fig = plot_rates(tm,ts,tu,tbreaks,Qrate,rs,scalefactor,yunits,xlab,ylab,titlestring,relative_dates)
% plot rates 
% input:
% tm             == master epoch in decimal years
% ts             == slave epoch in decimal years
% tbreaks        == break points in decimal years
% Qrate          == rates
% rs             == uncertainty on rates
% yunit          == string containing dimensions of rates
% xlab           == label for X-axis
% ylab           == label for Y-axis
% titelstring    == title string
% relative_dates == flag to use relative dates

% 20180115 Kurt Feigl

%% make relative dates?
if relative_dates == 1   
    t0 = min([tm; ts]); % initial epoch
    tm = tm-t0;
    ts = ts-t0;
    tbreaks = tbreaks - t0;
    tu = tu - t0;
end

% midpoints of each interval
tmid = (tm+ts)/2;

% half interval
th = (ts-tm)/2.;

%% Plot rates
h_fig = figure; 
set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
set(h_fig,'DefaultTextInterpreter','tex');
% axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
hold on;
% ind23m = find(tm >= dyear(2016, 03, 14) & tm <=  dyear(2016, 03, 24));
% ind23s = find(ts >= dyear(2016, 03, 14) & ts <=  dyear(2016, 03, 24));
indother = 1:numel(tm);

pest1 = mean(Qrate/scalefactor);
psig1 = std(Qrate/scalefactor);

errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'rs',5);
errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'rs',5);

ymin = nanmin(Qrate-rs)/scalefactor;
ymax = nanmax(Qrate+rs)/scalefactor;

% draw mean in blue
%errorbar_plus2((max(tu)+min(tu))/2.,pest1,(max(tu)-min(tu))/2.,psig1,'ko-',5);
tdura = max(tu)-min(tu);
tmean = min(tu)+tdura/2.;
errorbar_plus2(tmean,pest1,tdura/2.,psig1,'ko-',5);

for i=1:numel(tbreaks)
    plot([tbreaks(i) tbreaks(i)],[ymin ymax],'k--');
end
hold off;
ylabel(ylab);

% title(strcat(titlestring...
%     ,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits))...
%     ,'Interpreter','none');
% 20180627 correct units
title(strcat(titlestring...
    ,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1),psig1(1),yunits))...
    ,'Interpreter','none');
set(gca, 'FontSize', 14);
if relative_dates == 1
    xlabel(xlab);
else
    xlabel('Date');
    [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
    tick_labels = datetime(axy, axm, axd);
    set(gca, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75); %char(tick_labels)
end

return
end

