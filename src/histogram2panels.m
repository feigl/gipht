function histogram2panels(x1, x2, label1, label2, labely, title_string, fmt)
% make histograms of two random variables
% 2025/02/03 Kurt Feigl

xmin=min([reshape(x1,numel(x1),1);reshape(x2,numel(x2),1)],[],'all','omitnan');
xmax=max([reshape(x1,numel(x1),1);reshape(x2,numel(x2),1)],[],'all','omitnan');
%figure
subplot(2,1,1);
histogram(x1);
title(title_string);
%[m1,s1,H1]=plot_and_label_histogram(x1);
set(gca,'XLim',[xmin,xmax]);
xlabel(labely);
legend(label1,'location','northeast');
subplot(2,1,2);
histogram(x2);
%[m2,s2,H2]=plot_and_label_histogram(x2);
set(gca,'XLim',[xmin,xmax]);
xlabel(labely);
legend(label2,'location','northeast');
return
end