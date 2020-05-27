function ScoreTable= vertex_table(EdgeList, Scores0, doZscores, verbose)
% given n x 2 list of edges, make n x n table

narginchk(2,4);
nargoutchk(1,1);

if nargin < 3 || exist('doZscores','var') == 0
    doZscores = 0;
end
if nargin < 4 || exist('verbose','var') == 0
    verbose = 0;
end

[nEdges,ncols] = size(EdgeList);
if ncols ~= 2
    error('need list with 2 columns\n');
end
if nEdges < 1
    nEdges
    error('need list with at least 1 edge\n');
end


[EdgeList1, ia, ic] = unique(EdgeList,'rows');
ndup = numel(find(abs(ia-ic)>0));

% handle case of duplicate edges -- UNTESTED
if ndup > 0
    for i=1:numel(ic)
        Scores0(ic) = nanmean(Scores0(ic));
    end
    EdgeList = EdgeList(ia,:)
    Scores0 = Scores0(ia)
    warning('Duplicate Edges\n');
end


% find Vertices
Vertices = unique([EdgeList(:,1),EdgeList(:,2)]);
nVertices = numel(Vertices);
if nVertices < 2
    nVertices
    error('Need at least 2 vertices.\n');
end

Scores1 = nan(nVertices,nVertices);


% decide sign
if nanmin(colvec(abs(Scores0))) < 0
    signFactor = -1.0;   % this will make matrix anti-symmetric
else
    signFactor = +1.0;   % will preserve symmetry
end

% build matrix of scores
for j=1:nVertices
    for i=1:nVertices
        ii = find(EdgeList(:,2) == Vertices(i));
        jj = find(EdgeList(:,1) == Vertices(j));
        k = intersect(jj,ii);
        if numel(k) == 1
            Scores1(j,i) = Scores0(k);
            Scores1(i,j) = Scores0(k)*signFactor;  % make anti-symmetric
        end
    end
end

if doZscores == 1
    fprintf(1,'Table of Zscores\n');
    Tvals = zscore(Scores1);
else
    fprintf(1,'Table of Scores\n');
    Tvals = Scores1;
end

% decide format
maxabs = nanmax(colvec(abs(Tvals)));
minabs = nanmin(colvec(abs(Tvals)));
if minabs > 1. && maxabs < 100.
   rformat = '%#+8.2f'; % format for floating point numbers
elseif minabs > 10. && maxabs < 10000.
   rformat = '%#+8.0f'; % format for floating point numbers
elseif minabs > 0.001 && maxabs < 1. 
   rformat = '%#+8.4f'; % format for floating point numbers
else
   rformat = '%#+8.1E'; % format for small or large numbers
end

if doZscores == 1
    rformat = '%#+8.2f'; % format for floating point numbers
end



% label rows and columns of table
colnames= {};
rownames={};
for j=1:nVertices
    colnames{end+1} = sprintf('M%8d',Vertices(j));
    rownames{end+1} = sprintf('S%8d',Vertices(j));
end

% make a table
ScoreTable = array2table(Tvals);
ScoreTable.Properties.VariableNames = colnames;
ScoreTable.Properties.RowNames = rownames;


if verbose == 1
    fprintf('%9s ',' ');
    
    for j=1:nVertices
        fprintf(1,'%8d ',Vertices(j));
    end
    
    fprintf(1,'MaxAbsR_\n');
    fprintf('%9s ',' ');
    for j=1:nVertices+1
        fprintf(1,'%8s ','--------'); %
    end
    fprintf(1,'\n');
    
    for j=1:nVertices  % loop over rows
        for i=1:nVertices % loop over columns
            if i == 1
                fprintf(1,'%8d: ',Vertices(j));  % label row
            end
            if isfinite(Tvals(j,i)) == 1
                valstr = sprintf(rformat,Tvals(j,i));
            else
                valstr = '     .  ';
            end
            
            fprintf(1,'%8s ',valstr);
        end
        % print maxabs of row
        fprintf(1,strcat(rformat,'\n'),nanmax(abs(Tvals(j,:))));
    end
    fprintf(1,'%8s ','MaxAbsC_:');
    for i=1:nVertices
        fprintf(1,sprintf('%s ',rformat),nanmax(abs(Tvals(:,i))));
    end
    fprintf(1,sprintf('%s %s',rformat,'\n'),nanmax(abs(colvec(Tvals))));
end


return

end




