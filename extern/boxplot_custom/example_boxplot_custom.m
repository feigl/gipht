clear

%Examples :
X=rand(10,5);
X2=rand(100,5);
X3=rand(50,5);
C{1}=X;C{2}=X2;C{3}=X3;
list_labels={'Jan','Feb','Mar','Apr','May'};
legends={'Sample 1','Sample 2','Sample 3'};

% Simple Boxplots
figure('color','w')
subplot(2,4,1)
boxplot_custom(X,list_labels)
set(gca,'fontsize',10)
title('Simple boxplot')

%Simple boxplots with custom colors
subplot(2,4,2)
boxplot_custom(X(:,1),X2(:,3),X2(:,5),'linecolor','r','fillcolor','b')
set(gca,'fontsize',10)
title('Simple boxplot with custom colors')

% Grouped Boxplots
subplot(2,4,3)
boxplot_custom(X,X2,list_labels,'mode',2,'list_legends',legends(1:2))
set(gca,'fontsize',10)
title('Pairs of boxplots')
subplot(2,4,4)
boxplot_custom(X,X2,X3,list_labels,'mode',3,'list_legends',legends)
set(gca,'fontsize',10)
title('A group of 3 boxplots')

% Horizontal Boxplots
subplot(2,4,5)
boxplot_custom(X,list_labels,'direction',2)
set(gca,'fontsize',10)
title('Horizontal boxplots')

% Alternate y axis for the last column
subplot(2,4,6)
X4=rand(100,5);
X4(:,6)=nansum(X4,2);
list_labels{6}='Jan-May';
boxplot_custom(X4,list_labels,'axis_yy',1)
set(gca,'fontsize',10)
yyaxis left 
ylabel('Months')
yyaxis right 
ylabel('Total January-May')
title('Alternate y-axis for the last column')

% Add the mean
subplot(2,4,7)
boxplot_custom(X,X2,X3,list_labels,'mode',3,'list_legends',legends,'add_mean',1,'direction',2)
set(gca,'fontsize',10)
title('A group of 3 boxplots, horizontal, with the mean added')

%Changing the outlier's definition (from 1.5 IQR to 0.5)
subplot(2,4,8)
boxplot_custom(X,list_labels,'outlier_multiplier',0.5)
set(gca,'fontsize',10)
title('More strict outlier definition')