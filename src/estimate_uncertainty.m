
function msig = estimate_uncertainty(crit68,p0,p1,cost0,cost1,acosts1,atemps,atrials,lb,ub,pnames,objfun,fitfun,xdata,ydata,ysigma)
%function msig = estimate_uncertainty(crit68,p0,p1,cost0,cost1,acosts1,atemps,atrials,lb,ub,pnames,objfun,fitfun,xdata,ydata,ysigma)
% estimate uncertainties from width of cost valley
% Kurt Feigl 2010-MAR-15

% initialize variables
nslices = 100;
%p1vals = zeros(numel(acosts1),numel(p1)); % initialize  this matrix
msig = NaN * ones(numel(p1),1);            % intialize sigmas
ymin = cost1  - 0.05*abs(cost0-cost1);
ymax = cost0  + 0.05*abs(cost0-cost1);
costdiffthreshl = 5*abs(cost0-cost1)/nslices;
costdiffthreshr = 5*abs(cost0-cost1)/nslices;
adjthresh = 1.0e-5;
nf=0;

% loop over parameters
for i=1:numel(p1);
    lsig(i) = lb(i);leftisset=0;   % left  1-sigma limit
    rsig(i) = ub(i);rigtisset=0;   % right 1-sigma limit
    % remove underscore to avoid problem with axis label
    tmpstr = pnames{i};pname_no_=strrep(tmpstr,'_','\_');
    
    db = (ub(i)-lb(i));
     if db > abs(p1(i)-p0(i)) & numel(strfind(pnames{i},'Offset')) == 0
        fprintf(1,'Varying parameter: %s\n',pnames{i});
        if isfinite(acosts1(1)) == 1
            p1vals = atrials(:,i)';
        end
        for j = 1:round(1.1*nslices)
            % perturb only 1 parameter at a time from final value
            dpl(j) = p1(i) - (j-1) * db/nslices;
            p5l = p1; p5l(i) = dpl(j);
            %costs5l(j) = feval(objfun,p5l,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
            costs5l(j) = feval(objfun,p5l,fitfun,xdata,ydata,ysigma);
            
            %fprintf(1,'%5d %f %f\n',j,p5l(i),costs5l(j));
            % find left side of 68 percent confidence interval
            if costs5l(j) >= crit68  & dpl(j) <= p1(i) & leftisset==0
                lsig(i) = p5l(i);
                fprintf(1,'FOUND %5d %f %f %f\n',j,lsig(i),costs5l(j),crit68);
                leftisset=1;
            end
        end
        
        for j = 1:round(1.1*nslices)
            % perturb only 1 parameter at a time from final value
            dpr(j) = p1(i) + (j-1) * db/nslices;
            p5r = p1; p5r(i) = dpr(j);
            %costs5r(j) = feval(objfun,p5r,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
            costs5r(j) = feval(objfun,p5r,fitfun,xdata,ydata,ysigma);
            %fprintf(1,'%5d %f %f\n',j,p5r(i),costs5r(j));
            % find right side of 68 percent confidence interval
            if costs5r(j) >= crit68 & dpr(j) >= p1(i) & rigtisset==0
                rsig(i) = p5r(i);
                fprintf(1,'FOUND %5d %f %f %f\n',j,rsig(i),costs5r(j),crit68);
                rigtisset=1;
            end
        end
    end
    
    %    % find half-width of 68 percent confidence interval
    %    maximum
    %    if abs(lsig(i)-p1(i)) > abs(rsig(i)-p1(i))
    %       msig(i) = abs(lsig(i)-p1(i));
    %    else
    %       msig(i) = abs(rsig(i)-p1(i));
    %    end
    %    minimum
    %    if abs(lsig(i)-p1(i))  abs(rsig(i)-p1(i))
    %       msig(i) = abs(lsig(i)-p1(i));
    %    else
    %       msig(i) = abs(rsig(i)-p1(i));
    %    end
    %    average
    msig(i) = (abs(lsig(i)-p1(i)) + abs(rsig(i)-p1(i)) )/2;   %
    %    hypoteneuse
    %msig(i) = sqrt( (lsig(i)-p1(i))^2 + (rsig(i)-p1(i))^2 );   %
    
    %if (db > adjthresh * p0(i))
    if db > abs(p1(i)-p0(i)) && numel(strfind(pnames{i},'Offset')) == 0
        nf=nf+1;h(nf)=figure;
        
        dp = [dpl dpr];
        costs5 = [costs5l costs5r];
        % WARNING: Making 3x3 panels messes up this plot!
        %ipanel = ipanel+1;
        %subplot(npanrows,npancols,ipanel);
        %axis([-Inf Inf ymin ymax]);
        plot(dpl,costs5l,'r.-','MarkerFaceColor','r'); hold on;
        plot(dpr,costs5r,'r.-','MarkerFaceColor','r');
        
        %axis([-Inf Inf -Inf Inf]);
        if isfinite(acosts1(1)) == 1
            plot(p1vals,acosts1,'k+');
        end
        
        plot ([lb(i) lb(i)],    [ymin ymax],   'k:' ,'Clip','off','LineWidth',2);
        plot ([ub(i) ub(i)],    [ymin ymax],   'k:' ,'Clip','off','LineWidth',2);
        plot ([p1(i) p1(i)],    [cost1 ymax],  'k-' ,'Clip','off','LineWidth',2);
        plot ([p0(i) p0(i)],    [cost0 ymax],  'k--','Clip','off','LineWidth',2);
        plot ([min(dp) max(dp)],[cost0  cost0],             'k--','Clip','off','LineWidth',1);  % cost of initial model
        %plot ([min(dp) max(dp)],[crit95 crit95],            'k--','Clip','off','LineWidth',1);  % 95 percent confidence threshold
        plot ([min(dp) max(dp)],[crit68 crit68],            'k--','Clip','off','LineWidth',2);
        plot ([p1(i)-msig(i) p1(i)+msig(i)], [cost1  cost1],'b-' ,'Clip','off','LineWidth',4); % 1-sigma error bar at base
        %title '95 and 68 percent confidence, upper and lower bounds, sigma'
        %text(max(dp)-0.05*(max(dp)-min(dp)),crit95,sprintf(' 95 percent'),'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0],'HorizontalAlignment','Right');
        text(max(dp)-0.05*(max(dp)-min(dp)),crit68,sprintf(' 68 percent'),'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0],'HorizontalAlignment','Right');
        text(min(dp)+0.05*(max(dp)-min(dp)),cost0, sprintf(' initial')   ,'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0]);
        text(min(dp)+0.05*(max(dp)-min(dp)),cost1, sprintf(' final')     ,'BackgroundColor',[1 1 1],'Margin',5,'EdgeColor',[0 0 0]);
        %axis([-Inf Inf min(costs5) max(costs5)]);
        
        % Fight with Matlab to print the same number of digits after the decimal point.
        %       xticklabels=get(gca,'XTickLabel');xtickvals=str2num(xticklabels);
        %       for k=1:numel(xtickvals)
        %          xticklabels2{k}=sprintf('%9.3f',xtickvals(k));
        %       end
        %       set(gca,'XTickLabel',xticklabels2);
        yticklabels=get(gca,'YTickLabel');ytickvals=str2num(yticklabels);
        for k=1:numel(ytickvals)
            yticklabels2{k}=sprintf('%10.4f',ytickvals(k));
        end
        set(gca,'YTickLabel',yticklabels2);
        
        %xlabel(pname_no_);ylabel('circular mean deviation (cycles/datum)');
        xlabel(pname_no_);ylabel('cost');
        
        printpdf(sprintf('param%03d_%s.pdf',i,fitfun));
    end
end
return

