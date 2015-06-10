function h1 = plotpairs4(tm,ts,Qdiff,Qdsig,Qdmod ...
, tfit,Qdfit,Qfsigl,Qfsigu ...
, xlab,ylab,titlestring,tbreaks)
% function h1 = plotpairs5(tm,ts,Qdiff,Qdsig,Qdmod,pairnames ...
% , tfit,Qdfit,Qfsig ...
% , xlab,ylab,titlestring,tbreaks)
% 
% inputs:
% all these are column vectors with length equal to the number of pairs
%   tm           master epochs in decimal years
%   ts           slave  epochs in decimal years
%   Qdiff        observed values of differential quantity (NOT rate)
%   Qdsig        uncertainty of observed volume values
%   Qdmod        modeled values of differential quantity
%   pairnames    names of pairs as cell
% Column vectors 
%  tfit time epochs
%  Qfit differential quantity as modeled by fit
%  Qfsig uncertainty on model fit
%
% 2014-06-27 Kurt Feigl
% 2014-07-15 add tbreaks



error(nargchk(12,13,nargin));

ndat = numel(Qdiff);
%tmid = (tm+ts)/2;
tu = unique([tm ts]);


h1 = figure('color','w'); hold on; set(h1,'DefaultTextInterpreter','None');    

% if exist('tfit','var') == 1 && exist('ppredm','var') == 1 && exist('pest','var') == 1
%     fprintf(1,'Using incoming model\n');
% end

% draw green, vertical lines at times when slope breaks
if exist('tbreaks','var') == 1
    for i=1:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[min(Qfsigl) max(Qfsigu)],'g--');
    end
end



% plot observed values of differential change in red
if numel(Qdfit) > 0
    %fprintf(1,'i,tm(i),ts(i),ppredm(i),ppreds(i),Qdsig(i)\n');
    for i=1:ndat
%       plot observed differential quantity so that master epoch falls on modeled (fit) line
        ppredm(i) = interp1(tfit,Qdfit,tm(i),'linear');
        ppreds(i) = ppredm(i) +  Qdiff(i);
%        fprintf(1,'%3d %#10.4f %#10.4f %#12.4g %#12.4g %#12.4g\n',i,tm(i),ts(i),ppredm(i),ppreds(i),Qdsig(i));
        plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,'ro-','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
        plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)   ,ppreds(i)+Qdsig(i)],'b-', 'LineWidth',1);                                                                % error bar
        plot(ts(i)                               ,ppreds(i)-Qdsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
        plot(ts(i)                               ,ppreds(i)+Qdsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
    end
end

%plot modeled values  in black
%disp('tfit');size(tfit)
if numel(Qdfit) > 0
    %fprintf(1,'i,  tfit,  Qfit\n');
%     for i=1:numel(tfit)
%         fprintf(1,'%3d %#10.4f %#12.4g\n',i,tfit(i),Qdfit(i));
%     end
    plot(tfit,Qdfit,      'k-', 'LineWidth',2); % modeled value
    plot(tfit,Qfsigl,'k-.','LineWidth',1); % lower envelope 
    plot(tfit,Qfsigu,'k-.','LineWidth',1); % upper envelope 
end



axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
%axis fill;
%axis auto

set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','Xcolor','k','Ycolor','k');
h=title (titlestring); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','HorizontalAlignment','Center');
h=xlabel(xlab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
h=ylabel(ylab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');


return


