function h1 = plotpairs5_bradys(tm,ts, sat, Qdiff,Qdsig,Qdmod ...
, tfit,Qdfit,Qfsigl,Qfsigu, logy, logt ...
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



%error(nargchk(12,13,nargin));

ndat = numel(Qdiff);
%tmid = (tm+ts)/2;
tu = unique([tm ts]);


h1 = figure('color','w'); hold on; set(h1,'DefaultTextInterpreter','None');    

% if exist('tfit','var') == 1 && exist('ppredm','var') == 1 && exist('pest','var') == 1
%     fprintf(1,'Using incoming model\n');
% end

% % draw green, vertical lines at times when slope breaks
% if exist('tbreaks','var') == 1
%     for i=1:numel(tbreaks)
%         plot([tbreaks(i) tbreaks(i)],[min(Qfsigl) max(Qfsigu)],'g--');
%     end
% end

[l, m1,n1] = unique(sat,'first');
[c1,d1] =sort(m1);
l = l(d1);

%color = zeros(ndat,1);
count = zeros(numel(l),1);
for i = 1:ndat
    sat_i = sat{i};
    switch(sat_i)
        case {'ERS'} 
            color{i} = 'go-';
            count(1) = count(1)+1;
            if count(1) == 1
                ers_1 = i;
            end
        case {'ALS'}
            color{i} = 'mo-';
            count(2) = count(2)+1;
            if count(2) == 1
                als_1 = i;
            end
        case {'ENV'}
            color{i} = 'co-';
            count(3) = count(3)+1;
            if count(3) == 1
                env_1 = i;
            end
        case {'TSX'}
            color{i} = 'ro-';
            count(4) = count(4)+1;
            if count(4) == 1
                tsx_1 = i;
            end
    end
end

first_i = sort([ers_1, als_1, env_1, tsx_1]);

[l,m1,n1] = unique(color,'first');
[c1,d1] =sort(m1);
l = l(d1);
for i = 1:numel(l)
    csym = l{i};
    switch(csym)
        case {'go-'} 
            cleg{i} = 'ERS';
        case {'mo-'}
            cleg{i} = 'ALS';
        case {'co-'}
            cleg{i} = 'ENV';
        case {'ro-'}
            cleg{i} = 'TSX';
    end
end

%color = {'ro-','co-','bo-','go-'};color{c}


% plot observed values of differential change in red
if numel(Qdfit) > 0
    %fprintf(1,'i,tm(i),ts(i),ppredm(i),ppreds(i),Qdsig(i)\n');
     for i = 1:numel(first_i)
%         ppredm(i) = interp1(tfit,Qdfit,tm(first_i(i)),'linear');
%         ppreds(i) = ppredm(i) +  Qdiff(first_i(i));
%         plot([tm(first_i(i)),ts(first_i(i))],[ppredm(i)            ,ppreds(i)]         ,color{first_i(i)},'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
%plotting with midpoint between epochs on model fit   
        ppredmid(i) = interp1(tfit,Qdfit,mean([tm(first_i(i)) ts(first_i(i))]),'linear');
        ppredDelta(i) = Qdiff(first_i(i))/2;
        ppredm(i) = ppredmid(i) - ppredDelta(i);
        ppreds(i) = ppredmid(i) + ppredDelta(i);
        plot([tm(first_i(i)),ts(first_i(i))],[ppredm(i)            ,ppreds(i)]         ,color{first_i(i)},'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
        
    end
    legend(cleg)
    for i=1:ndat
% %       plot observed differential quantity so that master epoch falls on modeled (fit) line
%         ppredm(i) = interp1(tfit,Qdfit,tm(i),'linear');
%         ppreds(i) = ppredm(i) +  Qdiff(i);
%        % nx = linspace(ppredm(i), noise(i), 100);
%       %  sin_noise_plot(tm(i), ts(i), 1, ppredm(i), scalefactor, 'cos')
%         plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,color{i},'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
%          plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)   ,ppreds(i)+Qdsig(i)],'b-', 'LineWidth',1);                                                                % error bar
%         plot(ts(i)                               ,ppreds(i)-Qdsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
%         plot(ts(i)                               ,ppreds(i)+Qdsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
        
%plotting with midpoint between epochs on model fit
        ppredmid(i) = interp1(tfit,Qdfit,mean([tm(i) ts(i)]),'linear');
        ppredDelta(i) = Qdiff(i)/2;
        ppredm(i) = ppredmid(i) - ppredDelta(i);
        ppreds(i) = ppredmid(i) + ppredDelta(i);
       % nx = linspace(ppredm(i), noise(i), 100);
      %  sin_noise_plot(tm(i), ts(i), 1, ppredm(i), scalefactor, 'cos')
        plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,color{i},'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
         plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)   ,ppreds(i)+Qdsig(i)],'b-', 'LineWidth',1);                                                                % error bar
        plot(ts(i)                               ,ppreds(i)-Qdsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
        plot(ts(i)                               ,ppreds(i)+Qdsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
        
    end
end

% draw green, vertical lines at times when slope breaks
if exist('tbreaks','var') == 1
    for i=1:numel(tbreaks)
        plot([tbreaks(i) tbreaks(i)],[min(Qfsigl) max(Qfsigu)],'g--');
    end
end

%plot modeled values  in black
tfit_uns = tfit(tfit >= tbreaks(4) & tfit <= tbreaks(5));
Qdfit_uns = Qdfit(tfit >= tbreaks(4) & tfit <= tbreaks(5));
Qfsigl_uns = Qfsigl(tfit >= tbreaks(4) & tfit <= tbreaks(5));
Qfsigu_uns = Qfsigu(tfit >= tbreaks(4) & tfit <= tbreaks(5));
if numel(Qdfit) > 0
    plot(tfit,Qdfit,      'k-', 'LineWidth',2); % modeled value
    plot(tfit,Qfsigl,'k-.','LineWidth',1); % lower envelope 
   plot(tfit,Qfsigu,'k-.','LineWidth',1); % upper envelope 
end
scalevec = 1:10:132;
tfit_uns = tfit_uns(scalevec);
Qdfit_uns = Qdfit_uns(scalevec);
Qfsigl_uns = Qfsigl_uns(scalevec);
Qfsigu_uns = Qfsigu_uns(scalevec);
 plot(tfit_uns,Qdfit_uns,      'w--', 'LineWidth',3); % modeled value
    plot(tfit_uns,Qfsigl_uns,'k-.','LineWidth',1); % lower envelope 
   plot(tfit_uns,Qfsigu_uns,'k-.','LineWidth',1); % upper envelope 
   hold on 
plot(logt, -(logy./1000000).* 7.9148e+03, 'c-', 'LineWidth', 2);

%axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
%axis fill;
%axis auto
%axis ij
%set(gca,'Ydir','reverse')
set(gca,'Ydir','reverse','FontName','Helvetica','Fontsize',12,'FontWeight','bold','Xcolor','k','Ycolor','k');
h=title (titlestring); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','HorizontalAlignment','Center');
h=xlabel(xlab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
h=ylabel(ylab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');


return


