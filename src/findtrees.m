function [trees, Q, tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(t0,t1,dispflag,imast,mdate,islav,sdate)
% function [trees, Q, tepochs, iepochs, iuniqorbs, uniqdates] = findtrees(t0,t1,dispflag,imast,mdate,islav,sdate)
% given pair-wise list of epochs t0 and t1, return a list of species
% Previous Versions:
%    Kurt Feigl, CNRS
%   2007 MAY 09  Allow case of np = 2
%   2005 January version 1
%   2005 JUL 14 get unique dates and orbit numbers, too, please
%   2005 AUG 13 return:
%               ispecies    vector of indices to species for individual epochs 
%               jmast, jslav   vector of indices to species for pairs
%   [species] = findspecies(t0,t1)
%              t0 is a vector of epochs expressed in decimal years
%              t1 is a vector of epochs expressed in decimal years 
%            These 2 vectors must both have the same length
%   2007 DEC 07
%   handle case of np = 1
%   2009 DEC 12
%   do something with just 2 arguments
%
% Current version:
%   Specific to graph theory. 
%
% INPUTS: 
%   t0 -        vector of master epochs or chronologically first epochs in pairs (tm) 
%   t1 -        vector of slave epochs, or chronologically second epochs in pairs (ts)
%   disflag -   (optional) input for displaying information; 1 for yes (display) 0 for no;
%               default is 1 (yes)
%   imast -     (optional) indices of master epochs
%   mdate -     (optional) dates of master epochs
%   islav -     (optional) indices of slave epochs
%   sdate -     (optional) dates of slave epochs
% 
% OUTPUTS:
%   trees -     matrix containing indices of all epochs in distinct components
%               (components represented as rows, epochs represented as rows).  NaN is
%               assigned for elements where epoch corresponding to column is not in
%               component corresponding to row
%   Q     -     edge-vertex incidence matrix; rows represent pairs and columns
%               represent epochs.  
%               Q(i,j) = { -1  when jth epoch is master epoch of ith pair,
%                           1  when jth epoch is slave epoch of ith pair,
%                           0  otherwise}
%   tepochs -   vector of epochs
%   iepochs -   vector of indices corresponding to tepochs
%   iuniqorbs - vector of indices corresponding to unique orbits
%   uniqdates - vector of of unique dates
%
% Kurt L. Feigl, Elena C. Baluyut, UW-Madison 
% 2015-02-09



nargchk(2,7,nargin);
nargoutchk (4,6,nargout);

if nargin == 2
    dispflag = 1;
end% 

if dispflag == 1
    fprintf(1,'%s begins ...\n',mfilename);
end;

% initialize values to return
trees = NaN;
Q = NaN;
iepochs = NaN;
tepochs = NaN;

if length(t0) ~= length(t1) 
    error 'mismatch in length'
end

np = length(t0);
if dispflag == 1
    disp 'number of pairs'
    np
end

% if np < 2
%     error('Not enough pairs.')
%     return
% end

% Storing epochs
intconst = 10000;
i0 = round(intconst * t0);
i1 = round(intconst * t1);
[iepochs,iuniq,juniq] = unique([i0 i1]);
tepochs = iepochs/intconst;

% Unique epochs in years
tu=sort(unique([t0 t1]));

% Fill in some blanks
if nargin <= 3
    % find indices for master and slave
    imast = zeros(np,1);
    islav = zeros(np,1);
    for i=1:np
        imast(i) = find(abs(tu-t0(i)) < 1/366.0); % find master epoch index
        islav(i) = find(abs(tu-t1(i)) < 1/366.0); % find slave epoch index
    end
    
    % format strings for dates
    mdate = cell([np,1]);
    sdate = cell([np,1]);
    for i=1:np  
        [tyear,tmonth,tday] = dyear2yyyymmdd(t0(i)); mdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd'); % store master epoch (decimal years)
        [tyear,tmonth,tday] = dyear2yyyymmdd(t1(i)); sdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd'); % store slave epoch (decimal years)
    end    
end

% Find unique orbits and dates if desired
if nargout > 4
    ims = [imast islav];
    iuniqorbs = ims(iuniq); % find indices of unique orbits
    msdates = [mdate sdate]; char(msdates(:));
    uniqdates = msdates(iuniq); char(uniqdates(:)); 
end
    
% Set number of epoch indices
me = length(iepochs);

% Display values if desired
if dispflag == 1
    disp 'number of distinct epochs', me
    disp 'find differencing operator Q such that Q*epochs = pairs ...'
end

% Initialize edge-vertex incidence matrix Q
Q = zeros(np,me);
for p = 1:np
   for j = 1:me
      if i0(p) == iepochs(j)
         Q(p,j) = -1; % assign -1 if j_th epoch is the master epoch of pair p
      elseif i1(p) == iepochs(j) 
         Q(p,j) = 1; % assign +1 if j_th epoch is the slave epoch of pair p
      end
   end
end

% Find rank deficiency
rd = me - rank(Q);
rd0 = rd;
nspecies = rd;

if dispflag == 1
    disp 'Number of species'
    nspecies
end

% Store information on edge-vertex incidence matrix
nm = size(Q);
nrows = nm(1);

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
    N = null(Q);
elseif isequal(verLessThan('matlab', '15.0'), 0)
    warning('Check new version definition of null, may not be compatible with previous versions')
else
    % rational null space
    N = null(Q,'r');
end
  
[nrnull,ncnull] = size(N); % stores number of rows and columns of null matrix

if ncnull ~= nspecies | nrnull ~= me % size of null matrix should be equal to rank deficiency (and number of distinct components) and number of epochs
    error 'Counting problem'
else
   trees = NaN*ones(nspecies,me);  % set up matrix for storing components
    for j=1:nspecies
        k = 0;
        for i=1:me
           if abs(N(i,j)) > 0
                k = k+1;
                trees(j,k) = i; % store indices of epochs in jth component
            end
        end
    end
end

trees = sortrows(trees);

return

