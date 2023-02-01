function boxplot2(x1, x2, label1, label2, labely, title_string, fmt)
% make a box plot of two random variables
% https://www.mathworks.com/matlabcentral/answers/164929-how-to-do-a-boxplot-for-three-samples-of-different-sizes
% 2023/01/16 Kurt Feigl

% x1 = randn(600,1);
% x2 = randn(1440,1);
% x3 = randn(500, 1);
% x = [x1; x2; x3];
% g = [zeros(length(x1), 1); ones(length(x2), 1); 2*ones(length(x3), 1)];
% boxplot(x, g)

x = [x1; x2];
g = [1*ones(length(x1), 1); 2*ones(length(x2), 1)];
boxplot(x, g ...
    ,'BoxStyle','traditional' ...
    ,'Labels',{sprintf('%s (n = %d)',label1,numel(x1)), sprintf('%s (n= %d)',label2,numel(x2))} ...
    ,'Whisker',1.5 ...
    ,'Jitter',0.5 ...
    ,'Notch','on' ...
    ,'BoxStyle','outline');
hold on;

%% set up graphics in order of legend

set(findobj(gca,'tag','Outliers'),'Color','k','MarkerFaceColor','m','Marker','o','MarkerSize',3);
plot(nan,nan,'Color','k','MarkerFaceColor','m','Marker','o','MarkerSize',3,'LineStyle','none');
legstr{1}='Outlier';


% extrema whiskers 
% Multiplier for the maximum whisker length, specified as a positive
% numeric value. The maximum whisker length is the product of Whisker and
% the interquartile range.
% 
% boxplot draws points as outliers if they are greater than q3 + w × (q3 –
% q1) or less than q1 – w × (q3 – q1), where w is the multiplier Whisker,
% and q1 and q3 are the 25th and 75th percentiles of the sample data,
% respectively.
% 
% The default value for 'Whisker' corresponds to approximately +/–2.7σ and
% 99.3 percent coverage if the data are normally distributed. The plotted
% whisker extends to the adjacent value, which is the most extreme data
% value that is not an outlier.

% top whisker
L=findobj(gca,'tag','Upper Whisker');
set(L(1),'LineWidth',2,'LineStyle','-','Color','g');
set(L(2),'LineWidth',2,'LineStyle','-','Color','g');
L=findobj(gca,'tag','Upper Adjacent Value');
set(L,'LineWidth',2,'LineStyle','-','Color','g');
plot(nan,nan,'k','LineWidth',2,'LineStyle','-','Color','g');
legstr{end+1}='Whisker: 75th percentile + 1.5 * (interquartile range)';

% 75th percentile
plot(nan,nan,'b-','LineWidth',1');
legstr{end+1}='Top edge of blue box: 75th percentile.';
q75 = quantile(x1,0.75)
text(1.2,q75,sprintf(fmt,q75),'HorizontalAlignment','left' ,'VerticalAlignment','middle','FontSize',9); % 'BackgroundColor','w','Margin',1,'EdgeColor','b',
q75 = quantile(x2,0.75)
text(1.8,q75,sprintf(fmt,q75),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);

% top notch
%Two medians are significantly different at the 5% significance level if
%their intervals do not overlap. boxplot represents interval endpoints
%using the extremes of the notches or the centers of the triangular
%markers. The notch extremes correspond to q2 – 1.57(q3 – q1)/sqrt(n) and
%q2 + 1.57(q3 – q1)/sqrt(n), where q2 is the median (50th percentile), q1
%and q3 are the 25th and 75th percentiles, respectively, and n is the
%number of observations without any NaN values. If the sample size is
%small, the notches might extend beyond the end of the box.
set(findobj(gca,'tag','Notch'),'Color','b','LineWidth',1);
plot(nan,nan,'b<','LineWidth',1');
legstr{end+1}='Notch: 95% confidence interval for median = q_2 – 1.57(q_3 – q_1)/sqrt(n)';


% median
set(findobj(gca,'tag','Median'),'LineWidth',5);
plot(nan,nan,'r-','LineWidth',5');
legstr{end+1}='Central red mark: median. ';

% bottom notches
set(findobj(gca,'tag','Notch'),'Color','b','LineWidth',1);
plot(nan,nan,'b>','LineWidth',1');
legstr{end+1}='Notch: 95% confidence interval for median = q_2 + 1.57(q_3 – q_1)/sqrt(n)';

%
% 25th percentile
set(findobj(gca,'tag','Box'),'LineWidth',1);
plot(nan,nan,'b-','LineWidth',1');
legstr{end+1}='Bottom edge of blue box: 25th percentile.';
q25 = quantile(x1,0.25)
text(1.2,q25,sprintf(fmt,q25),'HorizontalAlignment','left' ,'VerticalAlignment','middle','FontSize',9); %'BackgroundColor','w','Margin',1,'EdgeColor','b',
q25 = quantile(x2,0.25)
text(1.8,q25,sprintf(fmt,q25),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);

 lower whisker
L=findobj(gca,'tag','Lower Whisker');
set(L(1),'LineWidth',2,'LineStyle','-','Color','g');
set(L(2),'LineWidth',2,'LineStyle','-','Color','g');
L=findobj(gca,'tag','Lower Adjacent Value');
set(L,'LineWidth',2,'LineStyle','-','Color','g');
plot(nan,nan,'k','LineWidth',2,'LineStyle','-','Color','g');
legstr{end+1}='Whisker: 25th percentile - 1.5 * (interquartile range)';


% outliers
set(findobj(gca,'tag','Outliers'),'Color','k','MarkerFaceColor','m','Marker','o','MarkerSize',3);
plot(nan,nan,'Color','k','MarkerFaceColor','m','Marker','o','MarkerSize',3,'LineStyle','none');
legstr{end+1}='Outlier';


% make the plot symmetric
axis([0.8,2.2,-Inf,Inf]);

% label
ylabel(labely);
title(title_string,'FontSize',12,'FontWeight','bold');
legend(legstr,'FontName','Times','FontSize',12,'Location','SouthOutside');

return
end