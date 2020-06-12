%% plot table from Comsol optimization
% 20200430 Kurt Feigl

% read the tables
% cd /Volumes/GoogleDrive/My Drive/comsol_trials4/LdM
%read_comsol_tables

close all
nf=0;


%% generate and read table
% collect variables from COMSOL runs
% read_comsol_tables_ldm
% CSV files contain only the long 32-character names for the trial variables, not the units or descriptions
%fileNameBIG = 'read_comsol_tables_BIG.csv'
% fileNameBIG = 'read_comsol_tables_ldm_BIG.csv'
% TBIG = readtable(fileNameBIG);
% MAT file contains units and descriptions
load('read_comsol_tables_ldm_BIG.mat')




% get initial start time
% ip = find(strcmp(Tparams.name,'t0')==1);
% t0 = Tparams.value(ip)
% y0a = t0/3600./24./365.25
% does not match Helene's value
y0 = 2007.0672 

% calculate dates in years
% first inflection point
y1 = y0 + years(seconds(table2array(TBIG(:,4)))); 
TBIG = [TBIG, table(y1,'VariableNames', {'y1'})];
TBIG.Properties.VariableUnits{'y1'} =  'year';
TBIG.Properties.VariableDescriptions{'y1'} = 'Date of first injection';

% start of second injection
y2 = y0 + years(seconds(table2array(TBIG(:,5)))); 
TBIG = [TBIG, table(y2,'VariableNames', {'y2'})];
TBIG.Properties.VariableUnits{'y2'} =  'year';
TBIG.Properties.VariableDescriptions{'y2'} = 'Date of second injection';

% sort
[ObjectiveVals,iSort] = sort(TBIG.Objective,'ascend');
TBIG = TBIG(iSort,:);


% get the names of the variables
varNames = TBIG.Properties.VariableNames;
[nTrials,nVariables] = size(TBIG);


    

% Q1(y1), Q2(y2), Q3(tend) are volumetric flow rates in m^3/s
% for this run tend = 2019.8
% final run will be tend = 2020.3


% print the optimal solution
Tbest=TBIG(1,:)

%% extract relevant columns
iCols = [1:3,6,11,12];
%iCols = [1:3,6]; % Q1, Q2, Q3, num

%% make scatter panels
nf=nf+1;figure;
%corrplot(TBIG,'varNames',TBIG.Properties.VariableNames);
[Rcorr,Pvalue,Handles] =corrplot(TBIG(:,iCols),'varNames',varNames(iCols));
% Try not to interpret underscores as Latex subscripts. Fails. 
% for i=1:numel(iCols)
%     for j=1:numel(iCols)
%         ax1 = get(Handles(i,j),'Parent');
%         ax1.XLabel.Interpreter='none';
%     end
% end
savefig(sprintf('%sFig%03d.fig',mfilename,nf));
print(gcf,'-dpdf',sprintf('%sFig%03d.pdf',mfilename,nf),'-r600','-fillpage','-painters');




%% make new table

%TBIG.Properties.VariableNames{i},TBIG.Properties.VariableDescriptions{i},TBIG.Properties.VariableUnits{i});

Toptimal = table(TBIG.Properties.VariableNames','VariableNames',{'name'});
Toptimal = [Toptimal, table(table2array(TBIG(1,:))','VariableNames',{'value'})];

% find statistics
% find critical value from quantiles
% obj68=quantile(ObjectiveVals,0.68);
% obj95=quantile(ObjectiveVals,0.95);

% For six degrees of freedom (Q1, Q2, Q3, t1, t2, num), then we should set the 68.3% confidence level at sqrt(7.04) = 2.65,
% (Chap 15.6 of Press, W. H., S. A. Teukolsky, B. P. Flannery, and W. T. Vetterling (1992), Numerical recipes in Fortran 77: volume 1, volume 1 of
% Fortran numerical recipes: the art of scientific computing, Cambridge university press.  052143064X).
 %  Chi-2 table on page table on page 692 
obj68 = sqrt(7.04) % 68 percent confidence from
obj95 = sqrt(12.8) % 68 percent confidence from

i68 = find(ObjectiveVals <= obj68);


x68mins = nan(nVariables,1);
x68maxs = nan(nVariables,1);
for i=iCols
    xvals = table2array(TBIG(:,i));
    x68min = nanmin(xvals(i68));
    x68max = nanmax(xvals(i68));   
    x68mins(i) = x68min;
    x68maxs(i) = x68max;
end
Toptimal = [Toptimal, table(x68mins,'VariableNames',{'x68min'})];
Toptimal = [Toptimal, table(x68maxs,'VariableNames',{'x68max'})];
Toptimal = [Toptimal, table(TBIG.Properties.VariableUnits','VariableNames',{'units'})];
Toptimal = [Toptimal, table(TBIG.Properties.VariableDescriptions','VariableNames',{'description'})];
writetable(Toptimal,sprintf('%s_Toptimal_%s.csv',mfilename,datestr(now,30)));



%% plot objective function for each parameter
% yvals = ObjectiveVals;
% ylab = 'Objective function [dimensionless]';
yvals = log10(ObjectiveVals);
ylab = 'log10(Objective)';

for i=iCols
    nf=nf+1;figure;
    hold on;
    xvals = table2array(TBIG(:,i));
      
    plot(xvals,yvals,'k+');
    plot(xvals(1),yvals(1),'kp','MarkerSize',16,'MarkerFaceColor','r');
        
    errorbar(xvals(1),yvals(1),[],[],xvals(1)-Toptimal.x68min(i),Toptimal.x68max(i)-xvals(1),'r');
    plot([Toptimal.x68min(i), Toptimal.x68max(i)], [log10(obj68), log10(obj68)],'k--');
    
    if ismember(i,[11,12])
        axis([2007, 2021, -1.7, 1.2]);
    elseif ismember(i,[6])
        axis([0, Inf, -1.7, 1.2]);
    else
        %axis([-Inf, Inf, -1.7, 1.2]);
        axis([0.9*nanmin([xvals; Toptimal.x68min(i)]), 1.1*nanmax([xvals; Toptimal.x68max(i)]), -1.7, 1.2]);
    end
    
    xlabel(sprintf('%s [%s]',TBIG.Properties.VariableNames{i}, TBIG.Properties.VariableUnits{i}));
    ylabel(ylab,'Interpreter','none');
    if i==6
        fmt = '%s %s %12.2E (%12.2E, %12.2E) [%s]';
    else
        fmt = '%s %s %12.4f (%12.4f, %12.4f) [%s]';
    end
    title(sprintf(fmt...
        ,TBIG.Properties.VariableNames{i} ...
        ,TBIG.Properties.VariableDescriptions{i} ...
        ,xvals(1)...
        ,Toptimal.x68min(i)...
        ,Toptimal.x68max(i)...
        ,TBIG.Properties.VariableUnits{i}));
    savefig(sprintf('%sFig%03d.fig',mfilename,nf));
    print(gcf,'-dpdf',sprintf('%sFig%03d.pdf',mfilename,nf),'-r600','-fillpage','-painters');
end
   

%% plot trade-off for each parameter
iCols = [1:3,6]; % Q1, Q2, Q3, num only - do not use absolute dates
for i=iCols
    for j=iCols
        if j > i
            nf=nf+1;figure;
            hold on;
            xvals = table2array(TBIG(:,i));
            yvals = table2array(TBIG(:,j));
             
            % build interpolation function
            %Fi = scatteredInterpolant(colvec(xvals),colvec(yvals),colvec(ObjectiveVals));
            Fi = scatteredInterpolant(colvec(xvals),colvec(yvals),colvec(log10(ObjectiveVals)+2)...
                ,'nearest','nearest');
            
            % build grid
            xvec = linspace(nanmin(xvals),nanmax(xvals),100);
            yvec = linspace(nanmin(yvals),nanmax(yvals),100);
            [XGRD,YGRD] = meshgrid(xvec,yvec);
            
            % interpolate values onto grid
            ZGRD = Fi(XGRD,YGRD);
%             % draw grid as image
%             % set the transparency for values with NaN
%             ALPHA=ones(size(ZGRD));
%             ALPHA(isnan(ZGRD))=0.0;
%             imagesc(xvec,yvec,ZGRD,'AlphaData',ALPHA);
            
            % draw contours
%            contour(XGRD,YGRD,ZGRD,2+[log10(obj68),log10(obj68)],'k-','LineWidth',1);
            %contour(XGRD,YGRD,ZGRD,2+[log10(obj95),log10(obj95)],'k--','LineWidth',1);
            
            %             % plot trial values
            %             plot(xvals,yvals,'k+');
            %             plot(xvals(i68),yvals(i68),'ko');
            
            % draw color coded circles, filled for significant ones
            pointSize = 128;
            scatter(xvals,yvals,pointSize,log10(ObjectiveVals)+2);
            scatter(xvals(i68),yvals(i68),pointSize,log10(ObjectiveVals(i68))+2,'filled');
            axis xy;
            axis tight;
            colormap('jet');
            
            
            % plot optimal value with error bars
            errorbar(xvals(1),yvals(1)...
                ,yvals(1)-Toptimal.x68min(j),Toptimal.x68max(j)-yvals(1)...
                ,xvals(1)-Toptimal.x68min(i),Toptimal.x68max(i)-xvals(1)...
                ,'kh','MarkerSize',16,'MarkerFaceColor','k','LineWidth',2);

            
            xlab = sprintf('%s (%s) [%s]',TBIG.Properties.VariableNames{i},TBIG.Properties.VariableDescriptions{i},TBIG.Properties.VariableUnits{i});
            ylab = sprintf('%s (%s) [%s]',TBIG.Properties.VariableNames{j},TBIG.Properties.VariableDescriptions{j},TBIG.Properties.VariableUnits{j});
            xlabel(xlab,'Interpreter','none');
            ylabel(ylab,'Interpreter','none');
            
            colorbar;
            title('log10(Objective) + 2');
            savefig(sprintf('%sFig%03d.fig',mfilename,nf));
            print(gcf,'-dpdf',sprintf('%sFig%03d.pdf',mfilename,nf),'-r600','-fillpage','-painters');
        end
    end
end
 