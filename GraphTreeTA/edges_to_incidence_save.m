function [Qiev, Trees] = edges_to_incidence(edges)
% function [Qiev, trees] = edges_to_incidence(edges)
% inputs:
%    edges == list of pairs of indices to vertices, one row per pair
% outputs: 
%     Qiev  == edge-vertex incidence matrix
%     trees == list of trees, giving list of vertices in each tree, 
%              one row per tree
%
% simplified from find_trees by Kurt L. Feigl, Elena C. Reinisch, UW-Madison
% 
% Reference:
% Reinisch, E.C., Cardiff, M. & Feigl, K.L. Graph theory for analyzing
% pair-wise data: application to geophysical model parameters estimated
% from interferometric synthetic aperture radar data at Okmok volcano,
% Alaska. J Geod 91, 9?24 (2017).
% https://doi.org/10.1007/s00190-016-0934-5
%
% 201200326


narginchk(1, 1);
nargoutchk(1, 2);

% verbose ?
dispflag = 1;


% if dispflag == 1
%     fprintf(1,'%s begins ...\n',mfilename);
% end

% initialize values to return
Qiev = NaN;
ivertices = NaN;
tepochs = NaN;

[nedges, ncols] = size(edges);
if ncols ~= 2
    ncols
    error(sprintf('input matrix edges must have only 2 columns. It has %d\n',ncols));
end

if dispflag == 1
    disp 'number of pairs'
    nedges
end

% 
%% integer indices to vertices
i1 = edges(:,1);
i2 = edges(:,2);


[ivertices,iuniq,juniq] = unique([i1 i2]);
 
% find indices
ibegins = zeros(nedges,1);
ifinish = zeros(nedges,1);
for i=1:nedges
    ibegins(i) = find(ivertices == i1(i)); % index of starting vertex
    ifinish(i) = find(ivertices == i2(i)); % index of finishing vertex
end
    
% Set number of vertices
nvertices = length(ivertices);

% Display values if desired
if dispflag == 1
    disp 'number of distinct epochs', nvertices
    disp 'find differencing operator Q such that Q*epochs = pairs ...'
end

% Initialize edge-vertex incidence matrix Qiev
Qiev = zeros(nedges,nvertices);
for p = 1:nedges
   for j = 1:nvertices
      if i1(p) == ivertices(j)
         Qiev(p,j) = -1; % assign -1 if j_th epoch is the master epoch of pair p
      elseif i2(p) == ivertices(j) 
         Qiev(p,j) = 1; % assign +1 if j_th epoch is the slave epoch of pair p
      end
   end
end

Qiev

% Find rank deficiency
rd = nvertices - rank(Qiev);

% number of trees is equal to the rank defiency by Theorem 1 
ntrees = rd;

if dispflag == 1
    fprintf(1, 'Number of trees ntrees = %d\n',ntrees);
end

% Store information on edge-vertex incidence matrix
[nrows,mcols] = size(Qiev);


% NULL	Nullspace of a matrix
%  	   N = NULL(A) uses the pivoting LU factorization computed by
%  	   PLU and the resulting reduced row echelon form computed by REF to
%  	   find a matrix N whose columns are a basis for the nullspace of A.
%  	   The number of columns of N is the nullity of A.
%  	   If A has independent columns, then N is empty, N = [];
%  
%  	   (This supersedes the MATLAB function NULL(A) which computes a
%  	   basis for the nullspace of A with orthonormal columns and, for
%  	   badly conditioned problems, even a possibly different dimension.)

% if MATLAB version is too old or too recent, this function may not work as
% intended
if verLessThan('matlab', '7.10')
    warning('Using old version of null')
    N = null(Qiev);
elseif isequal(verLessThan('matlab', '15.0'), 0)
    warning('Check new version definition of null, may not be compatible with previous versions')
else
    % rational null space
    N = null(Qiev,'r');
end
% 20180716 should work
%N = null(Qiev)
% 20200327 Line above fails in Matlab release 2018b!
  
[nrnull,ncnull] = size(N); % stores number of rows and columns of null matrix

%% size of null matrix should be equal to rank deficiency (and number of distinct components) and number of epochs
% 20200327 use cell array

if ncnull ~= ntrees || nrnull ~= nvertices
    error 'Counting problem'
else
    for j=1:ntrees
        inonnull = find(abs(N(:,j)) > 0);
        % get row vector of indices of vertices
        ivertices1 = rowvec(ivertices(inonnull));
        % store indices of epochs in jth component
        Trees{j} = ivertices1;     
    end
end

% make one tree per row
Trees = Trees';

% print out
for i=1:ntrees
    fprintf(1,'Tree%03d contains vertices: ',i);
    fprintf(1,'%2d ',Trees{i});
    fprintf(1,'\n');
end


% sort trees
%Trees = sortrows(Trees);

return
end



