function plot_tradeoffs1
% plot trade-offs between parameters
load trialvals.mat
whos

% count the number of free parameters
iadj=find(abs(p0-p1)>0)
mfree = numel(iadj) % number of rows in plot
nprowsmax = 3; % plot as separate figures if greater than this

% draw a histogram
figure;
hist(acosts1,100);hold on;
plot([crit69 crit69], [0 numel(acosts1)/20],'r--','LineWidth',2');
xlabel('cost [cycles]');
ylabel('number of trials');

% find the estimates with cost less than critical
iok = find(acosts1 <= crit69);
ndat = numel(iok)
fprintf(1,'Number of points (%d) with cost values less than critical (%.4f)\n',ndat,crit69);
if ndat < 10
    warning('Too few subcritical points. Using all points.\n');
    iok = 1:numel(acosts1);
end
 
% % find the estimates with cost less than critical
% iok = find(acosts1 <= crit69);
% % calculate correlation coefficients
% CC = nan(mfree,mfree);
% 
% for i=1:mfree
%     for j=1:i
%         % pointers to index in parameter vector
%         ip = iadj(i); % plot on X axis
%         jp = iadj(j); % plot on Y axis
%         % calculate correlation coefficient
%         cc=corrcoef(trials(iok,ip),trials(iok,jp));
%         cc=cc(1,2);
%         CC(i,j) = cc;
%     end
% end



% npanels = mfree*(mfree-1)/2;
% kpanel = 0;
figure;
if mfree < nprowsmax
    FontSize=5;
else
    FontSize=12;    
end
set(gca,'FontSize',FontSize);
for i=1:mfree
    for j=1:i
        % pointers to index in parameter vector
        ip = iadj(i); % parameter to plot on X axis
        jp = iadj(j); % parameter to plot on Y axis
        %figure; hold on;
        %kpanel = kpanel+1;
        %kpanel = (i-1)*mfree + j;
        %kpanel = (j-1)*mfree + i;
        kpanel = mfree^2 - (j-1)*mfree - i + 1;
        
        xlab=sprintf('%s',strrep(char(pnames{ip}),'_',' '))
        ylab=sprintf('%s',strrep(char(pnames{jp}),'_',' '))
        
        if mfree < nprowsmax
            subplot(mfree,mfree,kpanel);
            set(gca,'FontName','Helvetica','FontWeight','Normal','FontSize',12);
        else
            figure
            set(gca,'FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
        end
        hold on;
        
         
        
        if i==j
            % panel is on diagonal
            hist(trials(iok,ip),50);
            if mfree < nprowsmax
                title(xlab,'FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
            else
                xlabel(xlab,'FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
                ylabel('Number of trials','FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
                printpdf(sprintf('HISTOGRAM_%s.pdf',char(pnames{ip})));
            end
        else
            fprintf(1,'i = %d j = %d kpanel = %d ip = %d jp = %d\n',i,j,kpanel,ip,jp);
            
            if mfree < nprowsmax
                if i == mfree || j == 1
                    axis on
                else
                    axis off
                end
            else
                xlabel(xlab,'FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
                ylabel(ylab,'FontName','Helvetica','FontWeight','Normal','FontSize',FontSize);
            end
            
            % get values of parameters to plot
            xp = trials(iok,ip);
            yp = trials(iok,jp);
            
            % plot all the trials that are better than the critical value
            plot(xp,yp,'.k');
            
            % calculate correlation coefficient
            cc=corrcoef(xp,yp);
            cc=cc(1,2);
            cstring = sprintf('r = %6.2f',cc);
            if abs(cc) > 0.4
                title(strcat('** ',cstring,' **'));
            else
                title(cstring);
            end
            fprintf(1,'%s\n',cstring);
            
            
            % plot the convex hull around the trials
            khull = convhull(xp,yp);
            x_trial_hull=xp(khull);
            y_trial_hull=yp(khull);
            plot(x_trial_hull,y_trial_hull,'r-');
            
            % find the optimal estimate
            iopt = find(acosts1 <= min(nanmin(acosts1)));
            x_opt = trials(iopt,ip);
            y_opt = trials(iopt,jp);
            %fprintf(1,'%12.4e %12.4e\n',x_opt,y_opt);
            plot(x_opt,y_opt,'ro','MarkerSize',4,'MarkerFaceColor','r');
            printpdf(sprintf('TRADEOFF_%s_%s.pdf',char(pnames{ip}),char(pnames{jp})));
        end
        
    end
end
if mfree < nprowsmax
    printpdf(sprintf('%s.pdf',mfilename));
end
return;

