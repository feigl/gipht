function yValsMod = invert_incidence_matrix(yValsObs,Qiev,itrees)
%function yValsMod = invert_incidence_matrix(yValsObs,Qiev,itrees)
%
% given yValsObs, estimate yValsMod using edge-vertex incidence matrix Qiev
%
% Kurt Feigl CNRS 
% 2005 January
% 2005 JUL 14 add orbit numbers and dates as option
%
% 2007 DEC replace Bp with Ddop
% 20160524 use datetime timetags
% 20200328 change name of function

verbose = 0;
if verbose == 1
fprintf(1,'%s begins ...\n',mfilename);
end

narginchk(3,3);
nargoutchk(1,1);

% set design matrix equal to edge-vertex incidence matrix
GG = Qiev;

% data vector
data = yValsObs;

% number of data
ndata = length(yValsObs);
if ndata > 1
    
    % number of trees
    [ntrees,ndummy] = size(itrees);
    
    % number of pairs is the number of edges in the graph
    % number of epochs is the number of vertices in the graph
    [nedges,nvertices] = size(Qiev);
    
    % add constraining equations as lines to design matrix GG
    for j=1:ntrees
        itree1 = itrees(j,:);
        k=isfinite(itree1);
        mk = length(find(k == 1));
        itree1 = itree1(1:mk);
        ndata = ndata+1;
        for i = 1:nvertices
            GG(ndata,i) = 0;
        end
        for i = 1:mk
            GG(ndata,itree1(i)) = 1;
            data(ndata) = 0;
        end
    end
    
    % Length of data vector, including constraints;
    md = size(data);
    % 'Dimensons of design matrix, including constraints
    nmdd2 = size(GG);
    % 'Rank defiency, including constraints; 
    rd = nvertices - rank(GG);
    
    if verbose == 1
        fprintf(1,'Length of data vector, including constraints %d\n',md);
        fprintf(1,'Dimensons of design matrix, including constraints nmdd2 = %d %d\n',nmdd2);
        fprintf(1,'Rank defiency, including constraints %d\n',rd);
        fprintf(1,'Begin least squares adjustment...\n');
    end
    
    if rd <= 1
        if ndata > 2
            yValsMod = inv(GG' * GG) * (GG' * data);
        else
            yValsMod(1) = 0;  % arbitrarily make first the origin
            yValsMod(2) = data(1);
        end
    else
        if ndata > 2
            warning('Rank defiency persists. Using pseudoinverse!');
            yValsMod = pinv(GG' * GG) * (GG' * data);
        else
            yValsMod(1) = 0;  % arbitrarily make first the origin
            yValsMod(2) = data(1);
        end
    end
else
    yValsMod = nan;
    warning(sprintf('Not enough data: ndata = %d\n',ndata));
end
return
end


