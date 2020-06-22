function Hfigs = plot_objective_function(fileNameBIG,iobjcolumn,iCols,ilogPlot,iSteps)
%% plot a table of objective functions
% 20206009 Kurt Feigl

% read the tables
% cd /Volumes/GoogleDrive/My Drive/comsol_trials4/LdM
%read_comsol_tables

narginchk(1,5);

% objective function is in first column
if nargin < 2
    iobjcolumn = 1;
end



%close all
nf=0;

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
[nTrials,nVariables1] = size(TBIG);

if nargin < 3
    iCols = setdiff([1:nVariables1]',iobjcolumn);
end

% default is to take log10 of objective function
if nargin < 4
    ilogPlot = 1;
end

% default is to make all plots
if nargin < 5
    iSteps = 1:4;
end


% sort by Objective function
ObjectiveVals = table2array(TBIG(:,iobjcolumn));
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





%
% % Q1(y1), Q2(y2), Q3(tend) are volumetric flow rates in m^3/s
% % for this run tend = 2019.8
% % final run will be tend = 2020.3
%
%
% % print the optimal solution
%Toptimal=TBIG(1,:)
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

% find statistics
% find critical value from quantiles
% obj68=quantile(ObjectiveVals,0.68);
% obj95=quantile(ObjectiveVals,0.95);

% For six degrees of freedom (Q1, Q2, Q3, t1, t2, num) then we should set the 68.3% confidence level at sqrt(7.04) = 2.65,
% (Chap 15.6 of Press, W. H., S. A. Teukolsky, B. P. Flannery, and W. T. Vetterling (1992), Numerical recipes in Fortran 77: volume 1, volume 1 of
% Fortran numerical recipes: the art of scientific computing, Cambridge university press.  052143064X).
%  Chi-2 table on page table on page 692
% icdf('chi2',0.954,6) = 12.8191
% obj68 = sqrt(7.04) % 68 percent confidence from
% obj95 = sqrt(12.8) % 68 percent confidence from


% For eight degrees of freedom (t1, t2, P1, P3, tmodel0, logeta_ac_ratio, ac, ),
% then we should set the 68.3% confidence level at the
% icdf('chi2',0.683,8 ) = 9.3078
% multiply by optimal value of Objective function, assuming that it was based on a measurement uncertainty of 1.0 m?
% This assumes that the objective function is an L2 norm (i.e., variance with units of data squared)
%obj68 = icdf('chi2',0.683,8 ) * Tbest.Objective
ndof = numel(iCols)
obj68 = icdf('chi2',0.683,ndof) * ObjectiveVals(1)

%% find the indices of values below threshold
i68 = find(ObjectiveVals <= obj68);

x68mins = nan(nVariables1,1);
x68maxs = nan(nVariables1,1);
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
writetable(Toptimal,sprintf('%s_Toptimal_%s.csv',mfilename,datestr(now,30)));

%% make grid of scatter plots
if ismember(1,iSteps) == true && nUniqueTrials > 1
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
        %         zliml = 0.9*nanmin(zvals);
        %         zlimu = 1.1*nanmax(zvals);
        z68max = obj68;
    case 1
        zvals = log10(ObjectiveVals);
        zlab = 'log10(Objective)';
        %         zliml = nanmin(zvals) - 0.1
        %         zlimu = nanmax(zvals) + 0.1
        z68max = log10(obj68);
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
        if xvalr - xvall > 0
            
            % start figure
            nf=nf+1;figure;
            hold on;
            plot(xvals,zvals,'kx');
            plot(xvals(1),zvals(1),'kp','MarkerSize',16,'MarkerFaceColor','r');
            
            % plot error bars
            if (isfinite(xvall) == true) && (isfinite(xvalr) == true)
                if (abs(xvall-xvalr) > eps)   && (xvall > 0)  && (xvalr > 0)
                    errorbar(xvals(1),zvals(1),[],[],xvall,xvalr,'r');
                end
            end
            % plot horizontal line at 68 percent confidence
            plot([Toptimal.x68min(i), Toptimal.x68max(i)], [z68max, z68max],'r--');
            
            %% find convex hull
            x68vals = [xvals(i68); x68mins(i); x68maxs(i)];
            y68vals = [zvals(i68); z68max;     z68max];
            DT = delaunayTriangulation(x68vals,y68vals);
            khull = convexHull(DT);
            plot(x68vals(khull), y68vals(khull),'b-');
            
            
            
            %     if ilogPlot == 1
            %         axis([0.9*nanmin([xvals; Toptimal.x68min(i)]), 1.1*nanmax([xvals; Toptimal.x68max(i)]), yliml, ylimu]);
            %     end
            
            %xlabel(sprintf('%s [%s]',TBIG.Properties.VariableNames{i}, TBIG.Properties.VariableUnits{i}));
            xlabel(sprintf('%s',TBIG.Properties.VariableNames{i}));
            ylabel(zlab,'Interpreter','none');
            if max(abs(xvals)) > 1.e4
                fmt = '%s %12.2E (%12.2E, %12.2E)';
            else
                fmt = '%s %12.4f (%12.4f, %12.4f)';
            end
            title(sprintf(fmt...
                ,TBIG.Properties.VariableNames{i} ...
                ,xvals(1)...
                ,Toptimal.x68min(i)...
                ,Toptimal.x68max(i)));
            %         ,TBIG.Properties.VariableUnits{i}));
            %     ,TBIG.Properties.VariableDescriptions{i}));
            %savefig(sprintf('%sFig%03d.fig',mfilename,nf));
            Hfigs(nf) = gcf;
            print(gcf,'-dpdf',sprintf('%sFig%03d.pdf',mfilename,nf),'-r600','-fillpage','-painters');
        end
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
                
                % build interpolation function
                Fi = scatteredInterpolant(colvec(xvals),colvec(yvals),colvec(zvals)...
                    ,'nearest','nearest');
                
                
                % build grid
                xvec = linspace(nanmin(xvals),nanmax(xvals),100);
                yvec = linspace(nanmin(zvals),nanmax(yvals),100);
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
                scatter(xvals,yvals,pointSize,zvals);
                scatter(xvals(i68),yvals(i68),pointSize,zvals(i68),'filled');
                axis xy;
                axis tight;
                colormap('jet');
                
                
                % plot optimal value with error bars
                errorbar(xvals(1),zvals(1)...
                    ,zvals(1)-Toptimal.x68min(j),Toptimal.x68max(j)-zvals(1)...
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
                title(sprintf('Trade-off between parameter %d and parameter %d\n',i,j));
                
                
                %savefig(sprintf('%sTrade%03dvs%03d.fig',mfilename,ii,jj));
                print(gcf,'-dpdf',sprintf('%sTrade%03dvs%03d.pdf',mfilename,ii,jj),'-r600','-fillpage','-painters');
                Hfigs(nf) = gcf;
            end
        end
    end
end
return
end
