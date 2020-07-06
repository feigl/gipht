function Hfigs = plot_objective_function(fileNameBIG,ndof,iColObj,iCols,ilogPlot,iSteps)
%% plot a table of objective functions
% 20200624 Kurt Feigl
% ndof == number of degrees of freedom (ndata - mparams)
% read the tables
% cd /Volumes/GoogleDrive/My Drive/comsol_trials4/LdM
%read_comsol_tables

narginchk(4,6);

% reduced chi2 is in second column
if nargin < 3
    iColObj = 2;
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
varNames = TBIG.Properties.VariableNames
[nTrials,nVariables1] = size(TBIG)
% show first five rows
TBIG(1:5,:)


% number of free parameters
mparams = numel(iCols)

% default is to take log10 of objective function
if nargin < 5
    ilogPlot = 1;
end

% default is to make all plots
if nargin < 6
    iSteps = 1:4;
end
iSteps


%% make new table with initial values
Tinitial = table(TBIG.Properties.VariableNames','VariableNames',{'name'})
Tinitial = [Tinitial, table(table2array(TBIG(1,:))','VariableNames',{'value'})];
Tinitial
writetable(Tinitial,sprintf('%s_Tinitial.csv',strrep(fileNameBIG,'.csv','')));


%% collect values of objective function
ObjectiveVals = table2array(TBIG(:,iColObj));
objectiveTag = 'FtestOnChi2';

%% prune zeros
TBIG = TBIG((ObjectiveVals > 0.),:);

% keep first value
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





%% make new table with optimal values
%Toptimal=TBIG(1,:);
Toptimal = table(TBIG.Properties.VariableNames','VariableNames',{'name'});
Toptimal = [Toptimal, table(table2array(TBIG(1,:))','VariableNames',{'value'})];
Toptimal


%% find statistics
switch objectiveTag
    case 'Ftest'
        %% Use F test
        % Gordon, R. G., S. Stein, C. DeMets, and D. F. Argus (1987), Statistical tests for closure of plate motion
        % circuits, Geophysical Research Letters, 14, 587-590. https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/GL014i006p00587
        %         The F-ratio test is a standard statistical test used to compare variances of distributions. Whereas the chi-square test of plate circuit
        %         closure described above focuses on the consistency of the final best-fit estimates of relative Euler vectors with the a priori assigned
        %         errors, the F-ratio test of plate circuit closure described here focuses on the differences in the overall fit of two different models
        %         (Figure 1) to the data.
        %         The test is formulated using chi2s and is analogous to the test of additional terms widely used in curve fitting. If the values of
        %         sigma are consistently overestimated by a constant multiplicative factor, the value of F is unaffected, which is a potential
        %         advantage over the chi-square test.
        
        alpha = 1. - 0.682 % significance level
        % take ratio of variances initial:best
        variance_ratio = ObjectiveVal0/ObjectiveVals(1)
        % find critical value of ratio
        %variance_ratio_critical = ftest_critical(alpha,variance_ratio,ndof,ndof)  % will make a plot
        % number of degrees of freedom is the same for both samples
        ndof1 = ndof;
        ndof2 = ndof;
        % this is a one-sided test
        variance_ratio_critical = icdf('F',1-alpha,ndof1,ndof2)   % just return the critical value
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
    case 'FtestOnChi2'
        %% Use F test on Chi2
        alpha = 1. - 0.682 % significance level
        % number of degrees of freedom is the same for both samples
        ndof1 = ndof;
        ndof2 = ndof;
        % this is a one-sided test
        %        chiSquare1 = ObjectiveVals(1);% optimum estimate
        %         chiSquare2 = ObjectiveVal0;  % initial estimate, presumably larger
        %         obj68 = fcritical * ObjectiveVals(1)
        
        % 20200705 test requires chi-square, not reduced chi-square
        chiSquare1 = ObjectiveVals(1) * ndof;% optimum estimate
        chiSquare2 = ObjectiveVal0    * ndof;  % initial estimate, presumably larger
        [fcritical,H,test_string] = ftest_chi2(alpha,chiSquare1,chiSquare2,ndof1,ndof2)
        obj68 = fcritical * ObjectiveVals(1) * ndof
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
writetable(Toptimal,sprintf('%s_Toptimal.csv',strrep(fileNameBIG,'.csv','')));


%% plot objective function for each parameter
ilogPlot
switch ilogPlot
    case 0
        zvalsAll = ObjectiveVals;
        zlab = 'Objective function [dimensionless]';
        z68min = ObjectiveVals(1);
        z68max = obj68;
        zval0 = ObjectiveVal0;
    case 1
        zvalsAll = log10(ObjectiveVals);
        zlab = 'log10(Objective)';
        z68min = log10(ObjectiveVals(1));
        z68max = log10(obj68);
        zval0  = log10(ObjectiveVal0);
    otherwise
        error('unknown ilogPlot')
end

%% plot the trajectory of improvement
%% make grid of scatter plots
if ismember(0,iSteps) == true && nUniqueTrials > 1
    TBIG(1:10,:)
    nf=nf+1;figure;
    % sort by call count
    [nCallCount,isort] = sort(table2array(TBIG(:,2),'ascend'));
    plot(nCallCount,zvalsAll(isort),'r+-');
    xlabel('number of evaluations');
    ylabel(zlab);
    Hfigs(nf) = gcf;
    print(gcf,'-dpdf',sprintf('%s_Fig%02d.pdf',strrep(fileNameBIG,'.csv',''),nf),'-r600','-fillpage','-painters');
end

if ismember(1,iSteps) == true && nUniqueTrials > 1 && exist('corrplot') == 2 && numel(iCols) > 1
    nf=nf+1;
    try
        [Rcorr,Pvalue,Handles] =corrplot(TBIG(:,iCols),'varNames',varNames(iCols));
        Hfigs(nf) = gcf;
        print(gcf,'-dpdf',sprintf('%s_Fig%02d.pdf',strrep(fileNameBIG,'.csv',''),nf),'-r600','-fillpage','-painters');      
    catch ME
        ME
        warning('corrplot failed');
    end
end

%% make a histogram of objective function
if ismember(2,iSteps) == true
    nf=nf+1;
    histogram(zvalsAll);
    xlabel(zlab);
    ylabel('count');
    Hfigs(nf) = gcf;
    print(gcf,'-dpdf',sprintf('%s_Fig%02d.pdf',strrep(fileNameBIG,'.csv',''),nf),'-r600','-fillpage','-painters');
end


%% plot objective function for each parameter
ipanel=0;
if ismember(3,iSteps) == true
    for i=iCols
        fprintf(1,'\n\n***Plotting parameter %s\n',TBIG.Properties.VariableNames{i});       
        ipanel=ipanel+1;
        % values of estimated parameter
        xvals = table2array(TBIG(:,i));
        % prune to avoid error message from polyshape
        iprunez = find(isfinite(zvalsAll));
        xvals = xvals(iprunez);
        zvals = zvalsAll(iprunez);
        
        % start figure
        nf=nf+1;
        if mparams < 13
            doPanels = true;
            subplot(ceil(sqrt(mparams)),floor(sqrt(mparams)),ipanel);%,'align');
            plotTag = 'ALL';
            MarkerSize = 9;
            FontSize = 9;
            FontWeight = 'normal';
        else
            doPanels = false;
            figure;
            plotTag = TBIG.Properties.VariableNames{i};
            MarkerSize = 12;
            FontSize = 12;
            FontWeight = 'bold';
        end
        hold on;
        % plot all trials
        plot(xvals,zvals,'b+');
        
        %% find convex hull
        %       include all points
        DT = delaunayTriangulation(xvals,zvals);
        khull = convexHull(DT);
        % PolyHull = polyshape(xvals(khull), zvals(khull); % BAD, BAD, by default polyshape attempts to simplify.
        % It fails for some parameters and produces the following warning.
        %         Warning: Polyshape has duplicate vertices, intersections, or other inconsistencies that may produce
        %         inaccurate or unexpected results. Input data has been modified to create a well-defined polyshape.
        PolyHull = polyshape(xvals(khull), zvals(khull),'Simplify',false);
        % if all is well, then the next 2 lines should plot the same hull.
        plot(PolyHull,'FaceColor',[0,0,1],'FaceAlpha',0.1,'EdgeColor',[0,0,1],'EdgeAlpha',1.0);
        plot(xvals(khull), zvals(khull),'k--','LineWidth',2);
        
        
        %% Find uncertainties
        Toptimal.x68min(i) = nan;
        Toptimal.x68max(i) = nan;

        %% find intersection between threshold line and hull
        % make a horizontal line segment on the top of the region of 68 percent confidence
        LineSeg68 = [nanmin(xvals), z68max; nanmax(xvals), z68max];
        % plot horizontal line at 68 percent confidence
        plot(LineSeg68(:,1),LineSeg68(:,2),'r:','LineWidth',2);
        %% find the intersection of the hull and segment
        SegmentIntersection = intersect(PolyHull,LineSeg68);
              
        %% the 68 percent confidence interval is bounded by the left-most and right-most points in the intersection
        if ~isempty(SegmentIntersection) 
            plot(SegmentIntersection(:,1),SegmentIntersection(:,2),'-','Color',[1,0,1]);
            Toptimal.x68min(i) = nanmin(SegmentIntersection(:,1));
            Toptimal.x68max(i) = nanmax(SegmentIntersection(:,1));
        else
            fprintf(1,'WARNING: Hull and line segment do not intersect for parameter %s\n',TBIG.Properties.VariableNames{i});
            % make a rectangular box outlining the region of 68 percent confidence
            PolyBox68=polyshape([nanmin(xvals), nanmax(xvals), nanmax(xvals), nanmin(xvals)]...
                ,[z68min,        z68min,        z68max,        z68max],'Simplify',false);
            plot(PolyBox68,'FaceColor','r','FaceAlpha',0.1,'EdgeColor','r','EdgeAlpha',1.0);
            
            %% find the intersection of the hull and the box
            % [PolyIntersection,out] = intersect(PolyHull,PolyBox68)
            % GOTCHA: intersection does not include the edges. To avoid this issue, use 'KeepCollinearPoints' switch 
            PolyIntersection = intersect(PolyHull,PolyBox68,'KeepCollinearPoints',true);
            % the 68 percent confidence interval is bounded by the left-most and right-most points in the intersection
            if ~isempty(PolyIntersection.Vertices)
                plot(PolyIntersection,'FaceColor',[1,0,1],'FaceAlpha',0.1,'EdgeColor',[1,0,1],'EdgeAlpha',1.0);
                Toptimal.x68min(i) = nanmin(PolyIntersection.Vertices(:,1));
                Toptimal.x68max(i) = nanmax(PolyIntersection.Vertices(:,1));
            else
              fprintf(1,'WARNING: Hull and box do not intersect for parameter %s\n',TBIG.Properties.VariableNames{i});            
            end

            
            %% find the difference of the box and the hull
%             PolySubtraction = subtract(PolyBox68,PolyHull);
%             % the 68 percent confidence interval is bounded by the left-most and right-most points in the difference
%             if ~isempty(PolySubtraction.Vertices)
%                 plot(PolySubtraction,'FaceColor',[1,0,1],'FaceAlpha',0.1,'EdgeColor',[1,0,1],'EdgeAlpha',1.0);
%                 Toptimal.x68min(i) = nanmin(PolySubtraction.Vertices(:,1));
%                 Toptimal.x68max(i) = nanmax(PolySubtraction.Vertices(:,1));
%             else
%                 fprintf(1,'WARNING: Hull and box do not intersect for parameter %s\n',TBIG.Properties.VariableNames{i});
%             end
        end

        % half-width of error bar to left
        xvall = abs(xvals(1)-Toptimal.x68min(i));
        
        % half-width of error bar to right
        xvalr = abs(Toptimal.x68max(i)-xvals(1));
             
        % plot error bars
        if (isfinite(xvall) == true) && (isfinite(xvalr) == true)
            if (abs(xvall-xvalr) > eps)   && (xvall > 0)  && (xvalr > 0)
                errorbar(xvals(1),zvals(1),[],[],xvall,xvalr,'r');
            end
        end

        % plot initial value
        plot(Tinitial.value(i),zval0,'ko','MarkerSize',MarkerSize,'MarkerFaceColor','g');
        
        % plot optimal value
        plot(Toptimal.value(i),z68min,'kp','MarkerSize',MarkerSize,'MarkerFaceColor','r');
         
        % set font for numerical values on labels
        set(gca,'FontWeight',FontWeight,'FontSize',FontSize);
        
        
        xlim([nanmin(xvals),nanmax(xvals)]);
        %xlabel(sprintf('%s [%s]',TBIG.Properties.VariableNames{i}, TBIG.Properties.VariableUnits{i}));
        
        xlabel(sprintf('%s',TBIG.Properties.VariableNames{i}),'FontWeight',FontWeight,'FontSize',FontSize);
        ylabel(zlab,'Interpreter','none','FontWeight',FontWeight,'FontSize',FontSize);
        if Toptimal.value(i) > 1.E3
            fmt = '%s %10.2E\n(%10.2E, %10.2E)';
        else
            fmt = '%s %10.4f\n(%10.4f, %10.4f)';
        end
        title(sprintf(fmt...
            ,TBIG.Properties.VariableNames{i} ...
            ,Toptimal.value(i) ...
            ,Toptimal.x68min(i) ...
            ,Toptimal.x68max(i) ...
            )...
            ,'FontWeight',FontWeight,'FontSize',FontSize);
        %           ,objectiveTag ...
        %         ,TBIG.Properties.VariableUnits{i}));
        %     ,TBIG.Properties.VariableDescriptions{i}));
        %savefig(sprintf('%sFig%03d.fig',mfilename,nf));
        
        if ipanel == mparams || doPanels == false
            Hfigs(nf) = gcf;
            print(gcf,'-dpdf'...
                ,sprintf('%s_%s_Fig%03d%s.pdf',sprintf('%s_Fig%02d.csv',strrep(fileNameBIG,'.csv','')),objectiveTag,nf,plotTag)...
                ,'-r600','-fillpage','-painters');
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
                    ,sprintf('%s_%s_Trade%03dvs%03d.pdf',sprintf('%s_Fig%02d.csv',strrep(fileNameBIG,'.csv','')),objectiveTag,i,j)...
                    ,'-r600','-fillpage','-painters');
                Hfigs(nf) = gcf;
            end
        end
    end
end
return
end
