function h1 = plotpairs_gpsandinsar(mse, tm,ts,Qdiff,Qdsig,Qdmod ...
, tfit,Qdfit,Qfsigl,Qfsigu ...
, xlab,ylab,titlestring,tbreaks, ref_pts)
% function h1 = plotpairs5(tm,ts,Qdiff,Qdsig,Qdmod,pairnames ...
% , tfit,Qdfit,Qfsig ...
% , xlab,ylab,titlestring,tbreaks)
% 
% INPUTS:
% all these are column vectors with length equal to the number of pairs
%   tm -             master epochs in decimal years
%   ts -             slave  epochs in decimal years
%   Qdiff -          observed values of differential quantity (NOT rate)
%   Qdsig -          uncertainty of observed volume values
%   Qdmod -          modeled values of differential quantity
%   tfit  -          time epochs
%   Qfit  -          differential quantity as modeled by fit
%   Qfsigl/Qfsigu -  lower/upper uncertainty on model fit 
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


% Read arguments
error(nargchk(12,17,nargin)); % incorrect number of input arguments

if nargin < 15
    ref_pts = 'mid'; % default to plotting along midpoint of pairs
end

% Initialize
ndat = numel(Qdiff);
tmid = (tm+ts)/2;
tu = unique([tm ts]);

% Start figure
h1 = figure('color','w'); hold on; set(h1,'DefaultTextInterpreter','tex'); 

if isreal(Qfsigl) == 0 || isreal(Qfsigu) == 0
    Qfsigl 
    Qfsigu
end

% Draw green, vertical lines at times when slope breaks
% if exist('tbreaks','var') == 1
%     for i=1:numel(tbreaks)
%         plot([tbreaks(i) tbreaks(i)],[min(Qfsigl) max(Qfsigu)],'g--');
%     end
% end

% Plot observed values of differential change in red
if numel(Qdfit) > 0
    for i=1:ndat
        if ts(i) - tm(i) < 0.003
            marker_color = 'ro-';
        else
            marker_color = 'ro-';
        end
        if strcmp(ref_pts, 'mast') == 1 % plot with master on model fit
            ppredm(i) = interp1(tfit,Qdfit,tm(i),'linear'); % set displacement value at master along model fit
            ppreds(i) = ppredm(i) +  Qdiff(i);              % set displacement value at slave
            % plot displacement
            plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,marker_color,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
            % plot error bar on slave
            plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)   ,ppreds(i)+Qdsig(i)],'b-', 'LineWidth',1);                                                                % error bar on slave
            plot(ts(i)                               ,ppreds(i)-Qdsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(ts(i)                               ,ppreds(i)+Qdsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar          
            
        else % set midpoint to fall on modeled fit
              % plot error bars on both master and slave
                          ppredmid(i) = interp1(tfit,Qdfit,mean([tm(i) ts(i)]),'linear'); 
            ppredDelta(i) = Qdiff(i)/2; % set displacement 
            ppredm(i) = ppredmid(i) - ppredDelta(i);    % set displacement value at master
            ppreds(i) = ppredmid(i) + ppredDelta(i);    % set displacement value at slave

            if abs((ts(i) - tm(i))- 0.0027 ) < .001
            plot([tm(i),tm(i)],[ppredm(i)-Qdsig(i)/sqrt(mse)   ,ppredm(i)+Qdsig(i)/sqrt(mse)],'r-.', 'LineWidth',1);                                                % (error bar)/sqrt(2) on master
            plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)/sqrt(mse)   ,ppreds(i)+Qdsig(i)/sqrt(mse)],'r-.', 'LineWidth',1);                                                % (error bar)/sqrt(2) on slave
            plot(tm(i)                               ,ppredm(i)-Qdsig(i)/sqrt(mse) ,'rv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(tm(i)                               ,ppredm(i)+Qdsig(i)/sqrt(mse) ,'r^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
            plot(ts(i)                               ,ppreds(i)-Qdsig(i)/sqrt(mse) ,'rv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(ts(i)                               ,ppreds(i)+Qdsig(i)/sqrt(mse) ,'r^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
                % plot displacement
            plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         , 'Color', [0.5, 0.5, 0.5],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
  

            else
            plot([tm(i),tm(i)],[ppredm(i)-Qdsig(i)/sqrt(mse)   ,ppredm(i)+Qdsig(i)/sqrt(mse)],'b-.', 'LineWidth',1);                                                % (error bar)/sqrt(2) on master
            plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)/sqrt(mse)   ,ppreds(i)+Qdsig(i)/sqrt(mse)],'b-.', 'LineWidth',1);                                                % (error bar)/sqrt(2) on slave
            plot(tm(i)                               ,ppredm(i)-Qdsig(i)/sqrt(mse) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(tm(i)                               ,ppredm(i)+Qdsig(i)/sqrt(mse) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
            plot(ts(i)                               ,ppreds(i)-Qdsig(i)/sqrt(mse) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
            plot(ts(i)                               ,ppreds(i)+Qdsig(i)/sqrt(mse) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
              % plot displacement
            plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,'Color',[0.4940, 0.1840, 0.5560],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
  

            end
        end
    end
end

% Plot modeled values  in black
if numel(Qdfit) > 0
    plot(tfit,Qdfit,      'k-', 'LineWidth',2); % modeled value
    plot(tfit,Qfsigl,'k-.','LineWidth',1); % lower envelope 
    plot(tfit,Qfsigu,'k-.','LineWidth',1); % upper envelope 
end


% Adjust plotting format
axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
axis auto

set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','Xcolor','k','Ycolor','k');
h=title (titlestring); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','HorizontalAlignment','Center');
h=xlabel(xlab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
h=ylabel(ylab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');

curax = axis;

if exist('tbreaks','var') == 1
    for i=1:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[curax(3), curax(4)],'g--');
    end
end

return


