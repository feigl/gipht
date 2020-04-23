function itrees = find_trees_from_incidence(Qiev,verbose)
% function Qiev  = edges_to_incidence(edges)
% inputs:
%    Qiev     == edge-vertex incidence matrix
%    verbose  == 0 for quiet (default), 1 for verbose
% outputs: 
%     
%     irees == list of trees, giving list of indices to vertices in each tree, 
%              one row per tree
%     Trees == list of trees, giving list of vertices in each tree, 
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
% 2012003228


narginchk(0, 2);
nargoutchk(0, 1);

% verbose ?
if nargin < 2
    verbose = 0;
end


if verbose == 1
    fprintf(1,'Function %s begins ...\n',mfilename);
end


% Store information on edge-vertex incidence matrix
[nedges,nvertices] = size(Qiev);

if verbose == 1
    fprintf(1, 'Number of vertices = %d\n',nvertices);
    fprintf(1, 'Number of edges    = %d\n',nedges);
end


% Find rank deficiency
rd = nvertices - rank(Qiev);

% number of trees is equal to the rank defiency by Theorem 1 
ntrees = rd;

if verbose == 1
    fprintf(1, 'Number of trees    = %d\n',ntrees);
end




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
    error 'Counting problem.'
else
    %Trees=cell(ntrees,1);
    itrees = nan(ntrees,nvertices);
    ivertices = 1:nvertices;
    for j=1:ntrees
        inonnull = find(abs(N(:,j)) > 0);
        % get row vector of indices of vertices
        ivertices1 = rowvec(ivertices(inonnull));
        % store indices of epochs in jth component
        %Trees{j} = ivertices1; 
        n1 = numel(ivertices1);
        itrees(j,1:n1) = ivertices1;
        for k=n1+1:nvertices
           itrees(j,k)=nan;
        end
    end
end

%Trees = nan;

% % make one tree per row
% Trees = Trees';
% 
% % print out
% if verbose == 1   
%     for i=1:ntrees
%         fprintf(1,'Tree%03d contains vertices: ',i);
%         fprintf(1,'%d ',Trees{i});
%         fprintf(1,'\n');
%     end
% end
% 
% 
% % sort trees
% %Trees = sortrows(Trees);

if verbose == 1
    fprintf(1,'Function %s ended successfully.\n',mfilename);
end


return
end



