function h = plot_pairwise_correlations1(pvals,pnames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% make scatter plot
[ndat,mparam] = size(pvals);
mfree = 0;
XVALS = nan(ndat,mparam);
XVALS = nan(ndat,mparam);
for j = 1:mparam % columnes
    ptemp = pvals(:,j);
    ptemp = ptemp(find(isfinite(ptemp)==1));
    %ptemp = colvec(ptemp(find(abs(ptemp)>0)));
    %if numel(find(abs(ptemp)>0)) > 0
    if nanstd(ptemp)/abs(nanmean(ptemp)) > 1.0e-6
        mfree = mfree+1;
        pnames2{mfree} = pnames{j};
        fprintf(1,'%d %d %s\n',mparam,mfree,char(pnames{j}));
        for i=1:ndat
            fprintf(1,'%5d %12.4E\n',i,ptemp(i));
            XVALS(i,mfree) = ptemp(i)-nanmean(ptemp);
            YVALS(i,mfree) = ptemp(i)-nanmean(ptemp);
        end
        fprintf(1,'%5s %12.4E\n','mean ',nanmean(ptemp));
        fprintf(1,'%5s %12.4E\n','std  ',nanstd(ptemp));
    end
end
%pnames2

XVALS = XVALS(1:ndat,1:mfree);
YVALS = YVALS(1:ndat,1:mfree);


% GPLOTMATRIX  Scatter plot matrix with grouping variable.
%     GPLOTMATRIX(X,Y,G) creates a matrix of scatter plots of the columns of
%     X against the columns of Y, grouped by G.  If X is P-by-M and Y is
%     P-by-N, GPLOTMATRIX will produce a N-by-M matrix of axes.  If you omit
%     Y or specify it as [], the function graphs X vs. X.  G is a grouping
%     variable that determines the marker and color assigned to each point in
%     each matrix, and it can be a categorical variable, vector, string
%     matrix, or cell array of strings.  Alternatively G can be a cell array
%     of grouping variables (such as {G1 G2 G3}) to group the values in X by
%     each unique combination of grouping variable values.

% GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ,DOLEG,DISPOPT,XNAM,YNAM) specifies XNAM
%     and YNAM as the names of the X and Y variables.  Each must be a
%     character array or cell array of strings of the appropriate dimension.
%G = [1:mfree];
G = 1:ndat;
% G = cell(ndat,1);
% for i=1:ndat
%     G{i} = 'k.';
% end
doleg = 'off';
%doleg = [];
%dispopt = 'variable';
dispopt = 'hist';
%dispopt = 'on';
%dispopt = 'off';
%siz = 5;
siz = [];
% [h,ax,bigax] = gplotmatrix(XVALS,YVALS,G...
%     ,'k','.',siz,doleg,dispopt...
%     ,strrep(pnames2,'_',' ')...
%     ,strrep(pnames2,'_',' '));
[h,ax,bigax] = gplotmatrix(XVALS,[],[]...
    ,'k','.',siz,doleg,dispopt...
    ,strrep(pnames2,'_',' ')...
    ,strrep(pnames2,'_',' '));
for i=1:numel(ax)
    % tick labels
    set(ax(i),'FontName','Times','Fontsize',5,'FontWeight','normal');
    % labels on X axis
    h1=get(ax(i),'XLabel');
    set(h1,'FontName','Times','Fontsize',5,'FontWeight','normal');
    % labels on Y axis
    h1=get(ax(i),'YLabel');
    set(h1,'FontName','Times','Fontsize',5,'FontWeight','normal');
end

h = gcf;
end

