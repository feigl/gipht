function [trees, DD] = findtrees(t0,t1,score)
% given pair-wise list of epochs t0 and t1, return a list of trees
% Kurt Feigl, CNRS
% 2007 MAY 09  Allow case of np = 2
% 2005 January version 1
% 2005 JUL 14 get unique dates and orbit numbers, too, please
% 2005 AUG 13 return:
%               ispecies    vector of indices to trees for individual epochs 
%               jmast, jslav   vector of indices to trees for pairs
%
% [trees] = findspecies(t0,t1)
% 
%              t0 is a vector of epochs expressed in decimal years
%              t1 is a vector of epochs expressed in decimal years 
%            These 2 vectors must both have the same length
% 2007 DEC 07
% handle case of np = 1
% 2009 DEC 12
% do something with just 2 arguments

nargchk(3,3,nargin);
nargoutchk (2,2,nargout);


% initialize values to return
trees = NaN;
DD = NaN;

if length(t0) ~= length(t1) 
    error 'mismatch in length'
end

np = length(t0);

intconst = 10000;
i0 = round(intconst * t0);
i1 = round(intconst * t1);
[iepochs,iuniq,juniq] = unique([i0 i1]);
tepochs = iepochs/intconst;

% unique epochs in years
tu=sort(unique([t0 t1]));

% % fill in some blanks
% if nargin <= 3
%     %find indices for master and slave
%     imast = zeros(np,1);
%     islav = zeros(np,1);
%     for i=1:np
%         imast(i) = find(abs(tu-t0(i)) < 1/365.0);
%         islav(i) = find(abs(tu-t1(i)) < 1/365.0);
%     end
%     
%     % format strings for dates
%     mdate = cell([np,1]);
%     sdate = cell([np,1]);
%     for i=1:np
%         %     mdate{i} = sprintf('%5d',imast(i));
%         %     sdate{i} = sprintf('%5d',islav(i));
%         [tyear,tmonth,tday] = dyear2yyyymmdd(t0(i)); mdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd');
%         [tyear,tmonth,tday] = dyear2yyyymmdd(t1(i)); sdate{i} = datestr(datenum(tyear,tmonth,tday),'yyyymmdd');
%     end    
% end
% 
% 
% ims = [imast islav];
% msdates = [mdate sdate]; char(msdates(:));
% uniqdates = msdates(iuniq); char(uniqdates(:));
%     

me = length(iepochs);
    fprintf(1,'number of distinct epochs = %d\n', me)
    fprintf(1,'finding Incidence Matrix DD with %d rows and %d columns\n');
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

% rank defiency
rd = me - rank(DD);
ntrees = rd;
fprintf(1,'Number of trees is equal to the rank defiency of incidence matrix DD = %d\n',ntrees); 


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
if verLessThan('matlab', '7.10')
    warning('Using old version of null')
    N = null(DD);
else
    % rational null space
    N = null(DD,'r');
end
  
[nrnull,ncnull] = size(N);
if ncnull ~= ntrees || nrnull ~= me 
    error 'Counting problem'
else
   trees = NaN*ones(ntrees,me);  
    for j=1:ntrees
        k = 0;
        for i=1:me
           if abs(N(i,j)) > 0
                k = k+1;
                trees(j,k) = i;
            end
        end
    end
end

trees = sortrows(trees);
return

