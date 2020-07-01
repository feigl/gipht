function Hfigs = plot_objective_function(fileNameBIG,ndof,iColVarNormRes,iColChi2,iCols,ilogPlot,iSteps)
%% plot a table of objective functions
% 20200624 Kurt Feigl
% ndof == number of degrees of freedom (ndata - mparams)
% read the tables
% cd /Volumes/GoogleDrive/My Drive/comsol_trials4/LdM
%read_comsol_tables

narginchk(3,7);

% variance of sum of squares of normalized residuals is in first column
if nargin < 4
    iColVarNormRes = 1;
end

% chi2 is in second column
if nargin < 5
    iColChi2 = 2;
end


%close all
nf=0;
Hfigs=[];

% extract relevant columns
%iCols = [1:7,12,13,15];
%iCols = [6,7,14,15]; % debug on ratio
%iCols = [6,7]; % debug on ratio
%iCols = [1:13,15];



%% generate and read table
% collect variables from COMSOL runs
%read_comsol_tables_ldm
% CSV files contain only the long 32-character names for the trial variables, not the units or descriptions
%fileNameBIG = 'read_comsol_tables_BIG.csv'
% fileNameBIG = 'read_comsol_tables_ldm_BIG.csv'
TBIG = readtable(fileNameBIG);
% MAT file contains units and descriptions
%load('read_comsol_tables_ldm_2020MAY20gd_BIG.mat')

% get the names of the variables
varNames = TBIG.Properties.VariableNames;
[nTrials,nVariables1] = size(TBIG)

if nargin < 5
    iCols = setdiff([1:nVariables1]',iColVarNormRes);
end

% number of free parameters
mparams = numel(iCols)

% default is to take log10 of objective function
if nargin < 6
    ilogPlot = 1;
end

% default is to make all plots
if nargin < 7
    iSteps = 1:4;
end
iSteps

%% print the initial solution
Tinitial=TBIG(1,:)

%% collect values of objective function
if isfinite(iColVarNormRes) == true
    ObjectiveVals = table2array(TBIG(:,iColVarNormRes));
    objectiveTag = 'Ftest';
elseif isfinite(iColChi2) == true
    ObjectiveVals = table2array(TBIG(:,iColChi2));
    objectiveTag = 'chi2';
else
    error('Column for objective function is not specified.');
end
ObjectiveVal0 = ObjectiveVals(1)

%% sort by Objective function
[ObjectiveVals,iSort] = unique(ObjectiveVals,'sorted');
TBIG = TBIG(iSort,:);
[nUniqueTrials,nVariables2] = size(TBIG);
if nVariables1 == nVariables2
    nVariables = nVariables1;
else
    nVariables1
    nVariables2
    error('number of variables miscounted');
end



%% print the optimal solution
Toptimal=TBIG(1,:)
%
%
%
%
%% make new table
%
% %TBIG.Properties.VariableNames{i},TBIG.Properties.VariableDescriptions{i},TBIG.Properties.VariableUnits{i});
%
Toptimal = table(TBIG.Properties.VariableNames','VariableNames',{'name'});
Toptimal = [Toptimal, table(table2array(TBIG(1,:))','VariableNames',{'value'})];


%% find statistics
switch objectiveTag
    case 'Ftest'
        %% Use F test
        alpha = 1. - 0.682 % significance level
        % take ratio of variances initial:best
        variance_ratio = ObjectiveVal0/ObjectiveVals(1)
        % find critical value of ratio
        variance_ratio_critical = ftest_critical(alpha,variance_ratio,ndof,ndof)
        obj68 = variance_ratio_critical * ObjectiveVals(1)
    case 'chi2'
        % For six degrees of freedom (Q1, Q2, Q3, t1, t2, num) then we should set the 68.3% confidence level at sqrt(7.04) = 2.65,
        % (Chap 15.6 of Press, W. H., S. A. Teukolsky, B. P. Flannery, and W. T. Vetterling (1992), Numerical recipes in Fortran 77: volume 1, volume 1 of
        % Fortran numerical recipes: the art of scientific computing, Cambridge university press.  052143064X).
        %  Chi-2 table on page table on page 692
        % icdf('chi2',0.954,6) = 12.8191 (12.8 in table)

        %% 20200701 KF calculate ADDITIVE increase in chi2
        % Press et al. (1992) page 690, Theorem C: If a(j) is drawn from the universe of simulated data sets with actual parameters a0, then the
        % quantity DeltaChi2 = chi2(a(j)) - chi2(a0) is distributed as a chi2 distribution with M degrees of freedom. Here the chi2 are all evaluated using
        % the fixed (actual) values of the data set D0. This theorem makes the connection between particular values of DeltaChi2 and the fraction of the
        % probability distribution that they enclose as an M-dimensional region, i.e., the confidence level of the M-dimensional region.
        Delta_chi2 = icdf('chi2',0.683,mparams)
        obj68 =  Delta_chi2 + ObjectiveVals(1) 
    otherwise
        error('unknown objectiveTag');
end


%% find the indices of values below threshold
i68 = find(ObjectiveVals <= obj68);
if numel(i68) <= 1
    warning(sprintf('No solutions other than optimal found below threshold. Setting confidence interval to lowest part per 1000\n'));
    obj68 = quantile(ObjectiveVals,0.001)
    i68 = find(ObjectiveVals <= obj68);
    objectiveTag='Per1000'
end

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
%Toptimal = [Toptimal, table(TBIG.Properties.VariableUnits','VariableNames',{'units'})];
%Toptimal = [Toptimal, table(TBIG.Properties.VariableDescriptions','VariableNames',{'description'})];
Toptimal
writetable(Toptimal,sprintf('%s_Toptimal_%s.csv',mfilename,datestr(now,30)));

%% make grid of scatter plots
if ismember(1,iSteps) == true && nUniqueTrials > 1 && exist('corrplot') == 2
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
    %savefig(sprintf('%sFig%03d.fig',mfilename,nf));
    Hfigs(nf) = gcf;
    print(gcf,'-dpdf',sprintf('%sFig%03d.pdf',mfilename,nf),'-r600','-fillpage','-painters');
end

%% plot objective function for each parameter
switch ilogPlot
    case 0
        zvals = ObjectiveVals;
        zlab = 'Objective function [dimensionless]';
        z68min = ObjectiveVals(1);
        z68max = obj68;
        zval0 = ObjectiveVal0;
    case 1
        zvals = log10(ObjectiveVals);
        zlab = 'log10(Objective)';
        z68min = log10(ObjectiveVals(1));
        z68max = log10(obj68);
        zval0  = log10(ObjectiveVal0);
    otherwise
        error('unknown ilogPlot')
end


%% make a histogram of objective function
if ismember(2,iSteps) == true
    nf=nf+1;
    histogram(zvals);
    xlabel(zlab);
    ylabel('count');
    Hfigs(nf) = gcf;
    print(gcf,'-dpdf',sprintf('%sHistObj.pdf',mfilename),'-r600','-fillpage','-painters');
end


%% plot objective function for each parameter
if ismember(3,iSteps) == true
    for i=iCols
        % values of estimated parameter
        xvals = table2array(TBIG(:,i));
        % half-width of error bar to left
        xvall = xvals(1)-Toptimal.x68min(i);
        % half-width of error bar to right
        xvalr = Toptimal.x68max(i)-xvals(1);
        %         if xvalr - xvall <= eps
        %             warning(sprintf('No solutions other than optimal found below threshold for parameter %s . Setting confidence interval to [min,max]\n',Toptimal.name{i}));
        %             xvall = nanmin(xvals)
        %             xvalr = nanmax(xvals)
        %         end
        % start figure
        nf=nf+1;figure;
        hold on;
        % plot all trials
        plot(xvals,zvals,'kx');
        % plot optimal
        plot(xvals(1),zvals(1),'kp','MarkerSize',16,'MarkerFaceColor','r');
        % plot initial
        val0s = table2array(Tinitial);
        plot(val0s(i),zval0,'ko','MarkerSize',16,'MarkerFaceColor','g');
        
        % plot error bars
        if (isfinite(xvall) == true) && (isfinite(xvalr) == true)
            if (abs(xvall-xvalr) > eps)   && (xvall > 0)  && (xvalr > 0)
                errorbar(xvals(1),zvals(1),[],[],xvall,xvalr,'r');
            end
        end
        
        %% find convex hull
        % plot all points below threshold
        x68vals = [xvals(i68); x68mins(i); x68maxs(i)];
        y68vals = [zvals(i68); z68max;     z68max];
%         Plot only three points
%         x68vals = [nanmin(xvals(i68)); Toptimal.value(i); nanmax(xvals(i68)); ];
%         y68vals = [z68max;             z68min;            z68max];

        DT = delaunayTriangulation(x68vals,y68vals);
        khull = convexHull(DT);
        plot(x68vals(khull), y68vals(khull),'b-');
        
         % plot horizontal line at 68 percent confidence
        plot([nanmin(xvals), nanmax(xvals)], [z68max, z68max],'m:','LineWidth',2);
       
                
        %xlabel(sprintf('%s [%s]',TBIG.Properties.VariableNames{i}, TBIG.Properties.VariableUnits{i}));
        xlim([nanmin(xvals),nanmax(xvals)]);
        xlabel(sprintf('%s',TBIG.Properties.VariableNames{i}));
        ylabel(zlab,'Interpreter','none');
        if max(abs(xvals)) > 1.e4
            fmt = '%s %s %12.2E (%12.2E, %12.2E)';
        else
            fmt = '%s %s %12.4f (%12.4f, %12.4f)';
        end
        title(sprintf(fmt...
            ,objectiveTag ...
            ,TBIG.Properties.VariableNames{i} ...
            ,xvals(1)...
            ,Toptimal.x68min(i)...
            ,Toptimal.x68max(i)));
        %         ,TBIG.Properties.VariableUnits{i}));
        %     ,TBIG.Properties.VariableDescriptions{i}));
        %savefig(sprintf('%sFig%03d.fig',mfilename,nf));
        Hfigs(nf) = gcf;
        print(gcf,'-dpdf'...
            ,sprintf('%s_%s_Fig%03d%s.pdf',mfilename,objectiveTag,nf,TBIG.Properties.VariableNames{i})...
            ,'-r600','-fillpage','-painters');
    end
end



%% plot trade-off for each pair of parameters
Toptimal(iCols,:)


if ismember(4,iSteps) == true
    for ii=1:numel(iCols)
        for jj=1:numel(iCols)
            i = iCols(ii);
            j = iCols(jj);
            if j > i
                nf=nf+1;figure;
                hold on;
                xvals = table2array(TBIG(:,i));
                yvals = table2array(TBIG(:,j));
                
%                 % build interpolation function
%                 Fi = scatteredInterpolant(colvec(xvals),colvec(yvals),colvec(zvals)...
%                     ,'nearest','nearest');
%                 
%                 
%                 % build grid
%                 xvec = linspace(nanmin(xvals),nanmax(xvals),100);
%                 yvec = linspace(nanmin(zvals),nanmax(yvals),100);
%                 [XGRD,YGRD] = meshgrid(xvec,yvec);
%                 
%                 % interpolate values onto grid
%                 ZGRD = Fi(XGRD,YGRD);
%                 % draw grid as image
%                 % set the transparency for values with NaN
%                 ALPHA=ones(size(ZGRD));
%                 ALPHA(isnan(ZGRD))=0.0;
%                 imagesc(xvec,yvec,ZGRD,'AlphaData',ALPHA);
                
                % draw contours
                %            contour(XGRD,YGRD,ZGRD,2+[log10(obj68),log10(obj68)],'k-','LineWidth',1);
                %contour(XGRD,YGRD,ZGRD,2+[log10(obj95),log10(obj95)],'k--','LineWidth',1);
                
                %             % plot trial values
                %             plot(xvals,yvals,'k+');
                %             plot(xvals(i68),yvals(i68),'ko');
                
                % draw color coded circles, filled for significant ones
                pointSize = 128;
                scatter(xvals,yvals,pointSize,zvals);
                scatter(xvals(i68),yvals(i68),pointSize,zvals(i68),'filled');
                axis xy;
                axis tight;
                colormap('jet');
                
                
                % plot optimal value with error bars
%                 errorbar(xvals(1),zvals(1)...
%                     ,zvals(1)-Toptimal.x68min(j),Toptimal.x68max(j)-zvals(1)...
%                     ,xvals(1)-Toptimal.x68min(i),Toptimal.x68max(i)-xvals(1)...
%                     ,'kh','MarkerSize',16,'MarkerFaceColor','k','LineWidth',2);
%               corrected 20200701
                errorbar(xvals(1),yvals(1)...
                    ,yvals(1)-Toptimal.x68min(j),Toptimal.x68max(j)-yvals(1)...
                    ,xvals(1)-Toptimal.x68min(i),Toptimal.x68max(i)-xvals(1)...
                    ,'kh','MarkerSize',16,'MarkerFaceColor','k','LineWidth',2);
                
                if isfield(TBIG.Properties,'VariableDescriptions') == true && isfield(TBIG.Properties,'VariableUnits') == true
                    xlab = sprintf('%s (%s) [%s]',TBIG.Properties.VariableNames{i},TBIG.Properties.VariableDescriptions{i},TBIG.Properties.VariableUnits{i});
                    ylab = sprintf('%s (%s) [%s]',TBIG.Properties.VariableNames{j},TBIG.Properties.VariableDescriptions{j},TBIG.Properties.VariableUnits{j});
                else
                    xlab = sprintf('%s',TBIG.Properties.VariableNames{i});
                    ylab = sprintf('%s',TBIG.Properties.VariableNames{j});
                end
                xlabel(xlab,'Interpreter','none');
                ylabel(ylab,'Interpreter','none');
                
                Cbar=colorbar;
                Cbar.Label.String = zlab;
                title(sprintf('%s Trade-off between parameter %d and parameter %d\n',objectiveTag,i,j));
                
                
                %savefig(sprintf('%sTrade%03dvs%03d.fig',mfilename,ii,jj));
                print(gcf,'-dpdf'...
                    ,sprintf('%s_%s_Trade%03dvs%03d.pdf',mfilename,objectiveTag,i,j)...
                    ,'-r600','-fillpage','-painters');
                Hfigs(nf) = gcf;
            end
        end
    end
end
return
end
