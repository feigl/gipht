function [species, DD, tepochs, iepochs, iuniqorbs, uniqdates] = findspecies(t0,t1,dispflag,imast,mdate,islav,sdate)
% given pair-wise list of epochs t0 and t1, return a list of species
% Kurt Feigl, CNRS
% 2007 MAY 09  Allow case of np = 2
% 2005 January version 1
% 2005 JUL 14 get unique dates and orbit numbers, too, please
% 2005 AUG 13 return:
%               ispecies    vector of indices to species for individual epochs 
%               jmast, jslav   vector of indices to species for pairs
%
% [species] = findspecies(t0,t1)
% 
%              t0 is a vector of epochs expressed in decimal years
%              t1 is a vector of epochs expressed in decimal years 
%            These 2 vectors must both have the same length
% 2007 DEC 07
% handle case of np = 1
% 2009 DEC 12
% do something with just 2 arguments

nargchk(2,7,nargin);
nargoutchk (4,6,nargout);

if nargin == 2
    dispflag = 1;
end% 

if dispflag == 1
    fprintf(1,'%s begins ...\n',mfilename);
end;

% initialize values to return
species = NaN;
DD = NaN;
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

intconst = 10000;
i0 = round(intconst * t0);
i1 = round(intconst * t1);
[iepochs,iuniq,juniq] = unique([i0 i1]);
tepochs = iepochs/intconst;

% unique epochs in years
tu=sort(unique([t0 t1]));

% fill in some blanks
if nargin <= 3
    %find indices for master and slave
    imast = zeros(np,1);
    islav = zeros(np,1);
    for i=1:np
        imast(i) = find(abs(tu-t0(i)) < 1/365.0);
        islav(i) = find(abs(tu-t1(i)) < 1/365.0);
    end
    
    % format strings for dates
    mdate = cell([np,1]);
    sdate = cell([np,1]);
    for i=1:np
        %     mdate{i} = sprintf('%5d',imast(i));
        %     sdate{i} = sprintf('%5d',islav(i));
        [tyear,tmonth,tday] = dyear2yyyymmdd(t0(i)); mdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd');
        [tyear,tmonth,tday] = dyear2yyyymmdd(t1(i)); sdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd');
    end    
end


if nargout > 4
    ims = [imast islav];
    iuniqorbs = ims(iuniq);
    msdates = [mdate sdate]; char(msdates(:));
    uniqdates = msdates(iuniq); char(uniqdates(:));
end
    

me = length(iepochs);
if dispflag == 1
    disp 'number of distinct epochs', me
    disp 'find differencing operator DD such that DD*epochs = pairs ...'
end
DD = zeros(np,me);
for p = 1:np
   for j = 1:me
      if i0(p) == iepochs(j)
         DD(p,j) = -1;
      elseif i1(p) == iepochs(j) 
         DD(p,j) = 1;
      end
   end
end

% figure;
% spy (DD);
% title 'Differencing operator DD such that DD*epochs = pairs ';

%disp 'Rank defiency of DD'; 
rd = me - rank(DD);
rd0 = rd;
nspecies = rd;
if dispflag == 1
    disp 'Number of species'
    nspecies
end

nm = size(DD);
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
%if verLessThan('matlab', '7.9.0.529')
% if verLessThan('matlab', '7.10')
%     warning('Using old version of null')
%     N = null(DD);
% else
%     % rational null space
%     N = null(DD,'r');
% end

% Matlab version 8.5.0.197613 (R2015a)
% 
% null	Nullspace of a matrix
% N = null(A) uses the pivoting LU factorization computed by
% PLU and the resulting reduced row echelon form computed by REF to
% find a matrix N whose columns are a basis for the nullspace of A.
% The number of columns of N is the nullity of A.
% If A has independent columns, then N is empty, N = [];
% 
% (This supersedes the MATLAB function null(A) which computes a
%     basis for the nullspace of A with orthonormal columns and, for
%         badly conditioned problems, even a possibly different dimension.)
%         

N = null(DD);

  
[nrnull,ncnull] = size(N);
if ncnull ~= nspecies | nrnull ~= me 
    error 'Counting problem'
else
   species = NaN*ones(nspecies,me);  
    for j=1:nspecies
        k = 0;
        for i=1:me
           if abs(N(i,j)) > 0
                k = k+1;
                species(j,k) = i;
            end
        end
    end
end

species = sortrows(species);
return

