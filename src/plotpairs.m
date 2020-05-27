function plotpairs(tm,ts,Ydobs,Ydsig,Ydmod ...
, tfit,Ydfit,Yfsigl,Yfsigu ...
, xlab,ylab,titlestring,tbreaks, ref_pts, relative_dates ...
, xaxislims, yaxislims ...
, showfig)
% function plotpairs5(tm,ts,Ydiff,Ydsig,Ydmod,pairnames ...
% , tfit,Ydfit,Yfsig ...
% , xlab,ylab,titlestring,tbreaks)
% 
% INPUTS:
% all these are column vectors with length equal to the number of pairs
%   tm -             master epochs in decimal years
%   ts -             slave  epochs in decimal years
%   Ydobs -          observed values of differential quantity (NOT rate)
%   Ydsig -          uncertainty of observed volume values
%   Ydmod -          modeled values of differential quantity
%   tfit  -          time epochs
%   Yfit  -          differential quantity as modeled by fit
%   Yfsigl/Yfsigu -  lower/upper uncertainty on model fit 
%   xlab -           string input for x label
%   ylab -           sting input for y label
%   titlestring -    string input for title of figure
%   tbreaks -        vector of breaks in parameterization (reference epochs)
%   ref_pts -        optional string input specifying whether to plot with master ('mast') along the modeled fit or midpoint of the pair ('mid') along the modeled fit
%                    default is midpoint 'mid'
%
%  OUTPUT:
%   h1 -             function handle for figure
%
% 2014-06-27 Kurt Feigl
% Updates:
%   2014-07-15 add tbreaks
%   2015-07-27 add options for fit (ref_pts);  Elena C. Baluyut, UW-Madison
%   2018-01-15 add option for relative_dates
%   2020-04-23 make green line thicker
%   2020-04-24 do not return graphics handle


%% Read arguments
narginchk(12,18); % count input arguments

if nargin < 14
    ref_pts = 'mid'; % default to plotting along midpoint of pairs
end

if nargin < 15
    relative_dates = 0; % default to absolute dates
end



% Initialize
ndat = numel(Ydobs);

%% plot dates relative to initial epoch in decimal years
if relative_dates == 1
    t0 = min([tm; ts])
    tm = tm - t0;
    ts = ts - t0;
    tbreaks = tbreaks - t0;
    tfit = tfit - t0;
end

%% calculate other epochs
tmid = (tm+ts)/2;
tu = unique([tm ts]);

% set limits
if nargin < 16
    xaxislims = [min(tu)-0.1, max(tu)+0.1];
end

if nargin < 17
    yaxislims = [-Inf, +Inf];
end

if nargin < 18
    showfig = 1;
end


%% Start figure
if showfig == 1
    h1 = figure('color','w');
else
    % do not show
    h1 = figure('visible','off');
end

set(h1,'DefaultTextInterpreter','tex'); % or 'none'
hold on; 

if isreal(Yfsigl) == 0 || isreal(Yfsigu) == 0
    Yfsigl 
    Yfsigu
end

% Draw green, vertical lines at times when slope breaks
if exist('tbreaks','var') == 1
    for i=1:numel(tbreaks)
       % plot([tbreaks(i) tbreaks(i)],[min(Yfsigl) max(Yfsigu)],'g--','LineWidth',2);  % only goes to fit
       %plot([tbreaks(i) tbreaks(i)],[nanmin([Yfsigl;Ydobs-Ydsig]) nanmax([Yfsigu;Ydobs+Ydsig])],'g--','LineWidth',2);  % include ends of error bars
        plot([tbreaks(i) tbreaks(i)],[nanmin([Yfsigl;Ydobs]) nanmax([Yfsigu;Ydobs])],'g--','LineWidth',2);  % include ends of error bars
    end
end

% Plot observed values of differential change in red
if numel(Ydfit) > 0
    for i=1:ndat
        
        if strcmp(ref_pts, 'mast') == 1 % plot with master on model fit
            ppredm(i) = interp1(tfit,Ydfit,tm(i),'linear'); % set displacement value at master along model fit
            ppreds(i) = ppredm(i) +  Ydobs(i);              % set displacement value at slave
            % plot displacement
            plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,'ro-','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
            % plot error bar on slave
            plot([ts(i),ts(i)],[ppreds(i)-Ydsig(i)   ,ppreds(i)+Ydsig(i)],'b-', 'LineWidth',1);                                                                % error bar on slave
            plot(ts(i)                               ,ppreds(i)-Ydsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(ts(i)                               ,ppreds(i)+Ydsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar          
            
        else % set midpoint to fall on modeled fit
            ppredmid(i) = interp1(tfit,Ydfit,mean([tm(i) ts(i)]),'linear'); 
            ppredDelta(i) = Ydobs(i)/2; % set displacement 
            ppredm(i) = ppredmid(i) - ppredDelta(i);    % set displacement value at master
            ppreds(i) = ppredmid(i) + ppredDelta(i);    % set displacement value at slave
            % plot displacement
            plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,'ro-','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4,'LineWidth',2); % segment between master and slave        
            % plot error bars on both master and slave
            plot([tm(i),tm(i)],[ppredm(i)-Ydsig(i)/sqrt(2)   ,ppredm(i)+Ydsig(i)/sqrt(2)],'b-', 'LineWidth',1);                                                % (error bar)/sqrt(2) on master
            plot([ts(i),ts(i)],[ppreds(i)-Ydsig(i)/sqrt(2)   ,ppreds(i)+Ydsig(i)/sqrt(2)],'b-', 'LineWidth',1);                                                % (error bar)/sqrt(2) on slave
            plot(tm(i)                               ,ppredm(i)-Ydsig(i)/sqrt(2) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(tm(i)                               ,ppredm(i)+Ydsig(i)/sqrt(2) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
            plot(ts(i)                               ,ppreds(i)-Ydsig(i)/sqrt(2) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(ts(i)                               ,ppreds(i)+Ydsig(i)/sqrt(2) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
          
        end
    end
end

% Plot modeled values  in black
if numel(Ydfit) > 0
    plot(tfit,Ydfit, 'k-', 'LineWidth',2); % modeled value
    plot(tfit,Yfsigl,'k-.','LineWidth',1); % lower envelope 
    plot(tfit,Yfsigu,'k-.','LineWidth',1); % upper envelope 
end


% Adjust plotting format
%axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
%axis([min(tu)-0.1, max(tu)+0.1, -Inf, +Inf]);
axis([xaxislims(1), xaxislims(2), yaxislims(1), yaxislims(2)]);
%axis auto
axis manual
%axis tight
% set(gca,'XLimMode','auto');
% set(gca,'YLimMode','auto');
box on

set(gca,'FontName','Helvetica','Fontsize',10,'FontWeight','bold','Xcolor','k','Ycolor','k');
title(titlestring,'FontName','Helvetica','Fontsize',10,'FontWeight','bold','HorizontalAlignment','Center'...
    ,'interpreter','none');
% xlabel(xlab,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
% ylabel(ylab,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
fixlabels(xlab,'%.1f',ylab,'%.0f',10,10);

return


