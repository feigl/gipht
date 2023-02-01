%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate basic and special boxplots. Aside from the simple boxplots, you
% can create horizontal boxplots, alternate y axis, multiple bars, and
% ignore empty samples. You can also add the mean over the boxplot, and
% select the method for defining outliers.
%
% This program accepts vectors, matrixes and cells of data.
%
% Input :
% - data : Vector, Matrix with N columns, or cells with N elements. You
%          can also add N independants vectors/matrixes/cells and the
%          function will then create N boxplots
% - list_labels : cell of strings used to label the N boxplots (optional)
% - 'mode' : Integer. Create groups of boxplots. Default value is 1.
%            Current maximum of groups supported is 6. (optional)
% - 'list_legends' : cell of strings used to create the legend, if mode >1
%            (optional)
% - 'direction' : 1 : vertical (default), or 2 : horizontal (optional)
% - 'axis_yy' : Allows a second y axis on the right. The number
%		corresponds to the number of boxplots (form the right)
%	        which will use the second axis (optional)
% - 'ignore_empty' : If = 0 (default), empty vectors will generate
%		     empty boxplots. If = 1 empty vectors will be skipped
%            (optional)
% - 'linecolor' : color of the edges of the boxplots (rgb code or letter).
%                 Default is black (optional)
% - 'fillcolor' : color used to fill the boxplots (rgb code or letter). If
%                 the mode is greater than one, this input must be a cell
%                 (ex : {'r','g','m'} or {[1 0 0],[0 1 0],[0 0 1]}).
%                 Default is white when mode = 1. (optional)
% - 'outliercolor' : color for outlier (rgb code or letter). Default is
%                    blue (optional)
% - 'add_mean' : If = 1, the mean value of each sample will be added on the
%                boxplots. If = 0 (default), the mean doesn't appear 
%               (optional)
% - 'mean_marker' : If a marker for the mean is added, select the marker 
%                   type (ex : 'o','s','^', etc.). Default is a circle
%                   (optional)
% - 'mean_face_color' : If a marker for the mean is added, select the fill
%                   color (rgb code or letter). If the mode is greater than
%                   one, this input must be a cell (ex : {'r','g','m'} or
%                   {[1 0 0],[0 1 0],[0 0 1]}).Default is red when mode =1
%                   (optional). 
% - 'mean_edge_color' : If a marker for the mean is added, select the edge
%                       color (rgb code or letter). Default is black 
%                       (optional)
% - 'mean_size' : Size of the marker if the mean is displayed (optional).
%                 Default if 6.
% - 'outlier_method' : The classic method to identify outliers is to use
%                      the inter-quartile range (IQR = 75th percentile -
%                      25th percentile), but some also use standard
%                      deviations (Mean ± 3*Standard deviation). 'IQR'
%                      (default) will use the inter-quartile range, and
%                      'sigma' will use the standard deviation. The
%                      standard deviation is only recommended if the
%                      sample's distribution is normal.
%                      IQR :  Outlier > Q75+1.5*IQR and < Q25-1.5*IQR
%                      sigma : Outlier > mean+3*std and < mean-3*std
% - 'outlier_multiplier' : Specify the multiplier for the outlier
%                          definition. Default if 1.5 (if 'outlier_method'
%                          is 'IQR') or 3 (if 'outlier_method is 'sigma').
%
% Output :
% - H : handles of the boxplots
%
%
% %Examples :
% X=rand(10,5);
% X2=rand(100,5);
% X3=rand(50,5);
% C{1}=X;C{2}=X2;C{3}=X3;
% list_labels={'Jan','Feb','Mar','Apr','May'};
% legends={'Sample 1','Sample 2','Sample 3'};
% 
% % Simple Boxplots
% figure('color','w')
% subplot(2,4,1)
% boxplot_custom(X,list_labels)
% set(gca,'fontsize',10)
% title('Simple boxplot')
% 
% %Simple boxplots with custom colors
% subplot(2,4,2)
% boxplot_custom(X(:,1),X2(:,3),X2(:,5),'linecolor','r','fillcolor','b')
% set(gca,'fontsize',10)
% title('Simple boxplot with custom colors')
% 
% % Grouped Boxplots
% subplot(2,4,3)
% boxplot_custom(X,X2,list_labels,'mode',2,'list_legends',legends(1:2))
% set(gca,'fontsize',10)
% title('Pairs of boxplots')
% subplot(2,4,4)
% boxplot_custom(X,X2,X3,list_labels,'mode',3,'list_legends',legends)
% set(gca,'fontsize',10)
% title('A group of 3 boxplots')
% 
% % Horizontal Boxplots
% subplot(2,4,5)
% boxplot_custom(X,list_labels,'direction',2)
% set(gca,'fontsize',10)
% title('Horizontal boxplots')
% 
% % Alternate y axis for the last column
% subplot(2,4,6)
% X4=rand(100,5);
% X4(:,6)=nansum(X4,2);
% list_labels{6}='Jan-May';
% boxplot_custom(X4,list_labels,'axis_yy',1)
% set(gca,'fontsize',10)
% yyaxis left 
% ylabel('Months')
% yyaxis right 
% ylabel('Total January-May')
% title('Alternate y-axis for the last column')
% 
% % Add the mean
% subplot(2,4,7)
% boxplot_custom(X,X2,X3,list_labels,'mode',3,'list_legends',legends,'add_mean',1,'direction',2)
% set(gca,'fontsize',10)
% title('A group of 3 boxplots, horizontal, with the mean added')
% 
% %Changing the outlier's definition (from 1.5 IQR to 0.5)
% subplot(2,4,8)
% boxplot_custom(X,list_labels,'outlier_multiplier',0.5)
% set(gca,'fontsize',10)
% title('More strict outlier definition')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [H]=boxplot_custom(data,varargin)

hold all

%% Default parameters

% Number of grouped boxplots. 2 will display pairs of boxplots, 3 will
% display trios, etc. Maximum is 6, but it can be expanded eventually.
% Default is 1
mode=1;

% Labels for each boxplots/groups of boxplots
list_labels={};

% Legend, if mode > 1
list_legends={};

% Orientation of the boxplots. 1 : vertical (default), 2 : horizontal
direction=1;

% If greater than zero, a second y axis will be used on the right side. The
% number will indicate how many boxplots (from the right) will use this
% axis.
axis_yy=0;

% If this value is set to 1, empty boxplot will be skipped. If set to zero
% (default), an empty boxplot will be shown but the order will be
% preserved.
ignore_empty=0;

% Color of the edges of the boxplots
linecolor='k';

% Color inside the boxplots (Temporary empty. Will be defined later if no 
% custom value is provided)
fillcolor_temp=[];

% Color of the outliers
outliercolor='b';

% By default, the mean value of each boxplot doesn't appear. 
add_mean=0;
mean_marker=[];
mean_edge_color=[];
mean_face_color_temp=[];
mean_size=[];

% By default, the outlier definition is standard (inter-quartile range)
outlier_method='iqr';
outlier_multiplier=[];

% Width of boxplots. If the mode > 1, the width is reduced progressively
total_width=[0.5 0.25 0.17 0.12 0.1 0.08];

% Center of the boxplots
central_points{1}=0;
central_points{2}=[-0.175 0.175];
central_points{3}=[-0.215 0 0.215];
central_points{4}=[-0.24 -0.08 0.08 0.24];
central_points{5}=[-0.25 -0.125 0  0.125  0.25];
central_points{6}=[-0.26 -0.156 -0.052 0.052  0.156 0.26];

% Acquire the first vector/matrix of numerical data are convert it into
% elements of a "big" cell : Data. This cell will contain all numerical data used to
% build the boxplots.
if ismatrix(data) && isnumeric(data)
    Data(1).data=convert_cell(data);
end
count=2; %Increase the data counter.

if nargin>1
    skip=0;
    for k=1:length(varargin)
        if skip==1
            skip=0;
            continue
        end
        
        %% Acquire the different optional parameters
        if strcmpi(varargin{k},'list_labels')
            list_labels=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'mode')
            mode=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'list_legends')
            list_legends=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'direction')
            direction=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'axis_yy')
            axis_yy=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'ignore_empty')
            ignore_empty=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'linecolor')
            linecolor=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'fillcolor')
            fillcolor_temp=varargin{k+1};           
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'outliercolor')
            outliercolor=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'add_mean')
            add_mean=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'mean_marker')
            mean_marker=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'mean_edge_color')
            mean_edge_color=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'mean_face_color')
            mean_face_color_temp=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'mean_size')
            mean_size=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'outlier_method')
            outlier_method=varargin{k+1};
            skip=1;
            continue
        end
        
        if strcmpi(varargin{k},'outlier_multiplier')
            outlier_multiplier=varargin{k+1};
            skip=1;
            continue
        end
        
        % Acquire other data sets (vector/matrix) and add it to the "big"
        % data cell
        if ismatrix(varargin{k}) && isnumeric(varargin{k})
            Data(count).data=convert_cell(varargin{k});
            count=count+1; %Increase the data counter.
        end
        
        % Acquire other data sets (cell) and add it to the "big"
        % data cell
        if iscell(varargin{k}) && isnumeric(varargin{k}{1}) && ...
                ismatrix(varargin{k}{1})
            data_temp=varargin{k};
            [nl,nc]=size(data_temp);
            if nc>=nl
                Data(count).data=data_temp;
            else
                Data(count).data=data_temp';
            end
            count=count+1; %Increase the data counter.
        end
        
        % Liste of labels for each boxplots/groups of boxplots
        if iscell(varargin{k}) && ischar(varargin{k}{1})
            list_labels=varargin{k};
        end
        
    end
end

%% Select the width/centers/colors according to the selected mode.
width=total_width(mode);
centers=central_points{mode};

if isempty(fillcolor_temp)
    %Default filling colors
    if mode==1
        fillcolor{1}='w';
    else
        fillcolor{1,1}=[153 204 255]./256;
        fillcolor{2,1}=[255 153 153]./256;
        fillcolor{3,1}=[0 153 0]./256;
        fillcolor{4,1}=[255 128 0]./256;
        fillcolor{5,1}=[178 102 255]./256;
        fillcolor{6,1}=[102 255 255]./256;
    end
else
    %Custom filling colors
    if mode==1
        fillcolor{1}=fillcolor_temp;
    else
        if ~iscell(fillcolor_temp)
            error('The "fillcolor" input must be a cell if the mode > 1')
        elseif length(fillcolor_temp)<mode
            error('The length of the "fillcolor" input must be >= than the mode')
        end
        fillcolor=fillcolor_temp;
    end
end

%% If values are specified for mean_marker, mean size, mean_edge_color or mean_face_color but add_mean = 0, we guess that the user want to add the mean
 if (~isempty(mean_marker) || ~isempty(mean_size) || ~isempty(mean_edge_color) || ~isempty(mean_face_color_temp)) && add_mean==0
     add_mean=1;
 end
  
%% If the mean is added, select default parameters if they weren't
%  specified by the user
if add_mean==1
    if isempty(mean_marker)
        mean_marker='o';
    end
    
    if isempty(mean_edge_color)
        mean_edge_color='k';
    end
    
    if isempty(mean_size)
        mean_size=6;
    end
    
    if isempty(mean_face_color_temp)
        if mode==1
            mean_face_color{1,1}='r';
        else
            mean_face_color{1,1}=[0 102 204]./256;
            mean_face_color{2,1}=[255 0 0]./256;
            mean_face_color{3,1}=[0 102 0]./256;
            mean_face_color{4,1}=[153 76 0]./256;
            mean_face_color{5,1}=[76 0 153]./256;
            mean_face_color{6,1}=[0 153 153]./256;
        end
    else %Custom face colors for the mean
        if mode==1
            mean_face_color{1}=mean_face_color_temp;
        else
            if ~iscell(mean_face_color_temp)
                error('The "mean_face_color" input must be a cell if the mode > 1')
            elseif length(mean_face_color_temp)<mode
                error('The length of the "mean_face_color" input must be >= than the mode')
            end
            mean_face_color=mean_face_color_temp;
        end
    end
   
end


%% If the outlier method is IQR and the outlier_multiplier isn't specified,
% default is 1.5. If the outlier sigma method is selected, and the
% outlier_mulltiplier isn't specified, default is 3.
if strcmpi(outlier_method,'iqr')
    if isempty(outlier_multiplier)
        outlier_multiplier=1.5;
    end
elseif strcmpi(outlier_method,'sigma') || strcmpi(outlier_method,'std')
    if isempty(outlier_multiplier)
        outlier_multiplier=3;
    end   
else
    error('Unrecognized outlier method. Must be "iqr" or "sigma"')
end

    

% Identify potential problems
if mode>1
    if length(list_legends)~=mode && ~isempty(list_legends)
        error('The legend and the mode don''t have the same number of elements')
    end
    
    if length(Data)~=mode
        error('The number of datasets isn''t the same as the mode')
    end
    
    if ignore_empty==1
        warning('Ignoring empty boxplots might cause synchronism problems when the mode is larger than one')
    end
end

% If data is a cell, each element of the cell is incorporated in the "big"
% cell. If mode >1, each element of the "big" cell is a subgroup
if iscell(data) && length(data)==mode
    for k=1:mode
        if iscell(data{k})
            Data(k).data=data{k};
        elseif ismatrix(data{k}) && isnumeric(data{k})
            Data(k).data=convert_cell(data{k});
        end
    end
    % If mode >1 and the data cell has less element than the mode, data is the
    % first element of the "big" cell.
elseif iscell(data) && length(data)~=mode
    [nl,nc]=size(data);
    if nl<=nc
        Data(1).data=data;
    else
        Data(1).data=data';
    end
end

% If a second y axis is used, verify which boxplots will use this axis
if axis_yy>0
    if direction==2
        error('Warning! The  two y axis can''t be used currently with the horizontal orientation')
    end
    if mode>1
        nD=length(Data(1).data);
    else
        if length(Data)==1
            nD=length(Data(1).data);
        else
            nD=length(Data);
        end
    end
    limit_yy=nD-axis_yy+1;
    if limit_yy<2
        error('The number of boxplots which uses the second y axis is greater thant the total number of boxplots')
    end
end

% Initialize the counter
count=1;

%Drawing the boxplots, for each element of the "big" cell
for D=1:length(Data)
    if mode>1
        count=1; %Reinitialize the counter if there are multiple subgroups
    end
    
    %If there are two y axis, use the left one by default
    if axis_yy>0
        yyaxis left
        set(gca,'ycolor','k')
    end
    
    [~,nC]=size(Data(D).data);
    
    for C=1:nC
        
        x=Data(D).data{C};
        x(~isfinite(x))=nan;
        
        % Case empty boxplot
        if(sum(~isnan(x)))==0
            if ignore_empty==0
                count=count+1;
            else
                if mode>1
                    warning('Warning! Empty boxplots are ignored. This might cause synchronism problems if mode > 1')
                end
            end
            continue
        end
        
        %If there are two y axis, verify if the boxplot uses the second one
        if axis_yy>0
            if count==limit_yy
                yyaxis right
                set(gca,'ycolor','k')
            end
        end
        
        %Verify if the user has the "quantile" function available (Statistical Toolbox or other)
        try
            quantile(x,0.25);
        catch
           error('Statistical toolbox unavailable. We suggest you download the "Quantiles" function by David Ferreira (https://www.mathworks.com/matlabcentral/fileexchange/70279-quantiles) and place it in the current directory') 
        end
            
        %Compute the 25th, 50th and 75th percentile
        Q50=nanmedian(x);
        Q75=quantile(x,0.75);
        Q25=quantile(x,0.25);
        
        %Identify the outliers, according to the selected method
        if strcmpi(outlier_method,'iqr') %Outlier method: IQR
            %Interquartile range
            E25_75=Q75-Q25;
            %Verify if there are outliers
            outliers_top=x(x>(Q75+outlier_multiplier*E25_75));
            outliers_bottom=x(x<(Q25-outlier_multiplier*E25_75));
            %Max value which is not an outlier
            top=max(x(x<=(Q75+outlier_multiplier*E25_75)));
            %Min value which is not an outlier
            bottom=min(x(x>=(Q25-outlier_multiplier*E25_75)));
        else %Outlier method: Standard deviation
            outliers_top=x(x>(nanmean(x)+outlier_multiplier*nanstd(x)));
            outliers_bottom=x(x<(nanmean(x)-outlier_multiplier*nanstd(x)));
            %Max value which is not an outlier
            top=max(x(x<=(nanmean(x)+outlier_multiplier*nanstd(x))));
            %Min value which is not an outlier
            bottom=min(x(x>=(nanmean(x)-outlier_multiplier*nanstd(x))));
        end
        
        %If there are outliers, treat them apart from the other data
        [nl,nc]=size(x);
        if nl>nc
            outliers=[outliers_top ;outliers_bottom];
        else
            outliers=[outliers_top outliers_bottom];
        end
        
        %Compute the corners/centers of the boxplot
        if mode==1
            xl=count+centers(1)-width/2;
            xc=count+centers(1);
            xr=count+centers(1)+width/2;
            color=fillcolor{1};
        else
            xl=count+centers(D)-width/2;
            xc=count+centers(D);
            xr=count+centers(D)+width/2;
            color=fillcolor{D};
        end
        
        if direction==1 %Vertical case
            
            %Draw the boxplot
            H(D)=fill([xl,xr,xr,xl],[Q25 Q25 Q75 Q75],color,'EdgeColor',linecolor,'LineWidth',1,'LineStyle','-','marker','none');
            if top>Q75 %Case where max(x) > Q75
                line([xc xc],[Q75 top],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
                line([xl,xr],[top top],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            end
            if bottom<Q25 %Case where min(x) < Q25
                line([xc xc],[Q25 bottom],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
                line([xl,xr],[bottom bottom],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            end
            line([xl,xr],[Q50 Q50],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            
            %Draw the outliers
            plot((xc).*ones(length(outliers),1),outliers,'+','color',outliercolor,'MarkerSize',5)
            
            %If specified, draw the mean value
            if add_mean==1
                plot(xc,nanmean(x),'marker',mean_marker,'markeredgecolor',mean_edge_color,'markerfacecolor',mean_face_color{D},'markersize',mean_size)
            end
            
        elseif direction==2 % Horizontal case
            
            %Draw the boxplot
            H(D)=fill([Q25 Q25 Q75 Q75],[xl,xr,xr,xl],color,'EdgeColor',linecolor,'LineWidth',1,'LineStyle','-','marker','none');
            if top>Q75
                line([Q75 top],[xc xc],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
                line([top top],[xl,xr],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            end
            if bottom<Q25
                line([Q25 bottom],[xc xc],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
                line([bottom bottom],[xl,xr],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            end
            line([Q50 Q50],[xl,xr],'Color',linecolor,'LineWidth',1,'LineStyle','-','marker','none')
            
            %Draw the outliers
            plot(outliers,(xc).*ones(length(outliers),1),'+','color',outliercolor,'MarkerSize',5)
            
            %If specified, draw the mean value
            if add_mean==1
                plot(nanmean(x),xc,'marker',mean_marker,'markeredgecolor',mean_edge_color,'markerfacecolor',mean_face_color{D},'markersize',mean_size)
            end
        end
        
        %Increment the counter
        count=count+1;
        
        
    end
end

% Write the labels for each boxplots/groups of boxplots
if direction==1
    set(gca,'xtick',1:(count-1),'fontsize',20)
    if ~isempty(list_labels)
        set(gca,'xticklabel',list_labels)
    end
    xlim([0 count])
elseif direction==2
    set(gca,'ytick',1:(count-1),'fontsize',20)
    if ~isempty(list_labels)
        set(gca,'yticklabel',list_labels)
    end
    ylim([0 count])
end

% Write the legend, if mode >1
if mode>1 && ~isempty(list_legends)
    legend(H,list_legends,'location','best','autoupdate','off')
end

box on

end



%Function to convert vector/matrix data into a cell element
function [data_cell]=convert_cell(data)

[nl,nc]=size(data);

if nl==1
    data=data';
end

data_cell=cell(1,nc);

for C=1:nc
    data_cell{1,C}=data(:,C);
end

end