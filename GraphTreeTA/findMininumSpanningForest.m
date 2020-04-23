function [i_msf, i_rep, nCycles, QievMSF, TreesMSF, WeightsMSF] = findMininumSpanningForest(Qiev, Trees, Weights, scalefactor)
% function [i_msf, i_rep, Nc ] = plotMSF( DD, trees, rs, scalefactor)
% Using the incidence matrix Q, trees, epochs, and rates, returns
% the minimim spanning forest of the full data set using Prim's Algorithm with edges weighted by the
% uncertainty of the data.  In this graph, nodes represent epochs and edges
% represent pairs (interferograms).  
%
% INPUTS:
%   Qiev        - edge-vertex incidence matrix (nedges by nvertices)
%   Trees       - matrix of components in data set
%   Weights     - weighting for edges, e.g. uncertainties of pairwise data 
%   scalefactor - scaling value for data set
%
% OUTPUTS:
%   i_msf       - indices of pairs (edges) in minimum spanning forest - 
%   i_rep       - indices of pairs (edges) eliminated from minimum spanning forest
%   Nc          - vector the length of ntrees showing the number of cycles
%                 per each component
%
% Elena C. Baluyut, UW-Madison  2015-02-20
% 20200326  Kurt Feigl, update comments, change name of function

narginchk(1,4);
nargoutchk(1,6);

[nedges,nvertices] = size(Qiev)

% Find pairs
pairlst = pairlist(Qiev);
[PN] = pairlist_tree(pairlst, Trees, Weights/scalefactor);
[ntrees maxvertices] = size(Trees)
fprintf(1,'PN is matrix containing information of pairs\n');
fprintf(1,'Itree Ipair is it rsofpair\n')
 
PN

% Round weights to nearest 100th
rRSs = Weights/scalefactor*100;
rRSs = round(rRSs);
Wrs = rRSs/100;

%Find mimnimum spanning tree (Prim)
[ i_msf, i_rep, nCycles ] = MSF(ntrees, Trees, PN);
i_r = []; % initialize vector storing indices of epochs to be removed

for jj = 1:numel(i_rep)
    [i_n mdum] = find(PN(:,2) == i_rep(jj)); % find rows of PN corresponding to epochs to be removed
    i_r(end+1,:) = i_n; % store indices of epochs to be removed
end

% Define PN matrix for MSF
PN_msf = PN;
PN_msf(i_r, :) = []; % remove rows corresponding to epochs not in MSF
max_end = max([max(PN_msf(:,3)) max(PN_msf(:,4))]); % find highest index included in MSF

% Initialize and build connectivity graph of MSF
Cf = zeros(max_end, max_end); % set up connectivity matrix according to number of epochs in MSF
mn = max_end; 
[hh dum] = size(PN_msf); % get number of pairs in MSF
for h = 1:hh
    Cf(PN_msf(h,3), PN_msf(h,4)) = Weights(h)/scalefactor; % store uncertainty for each pair in MSF
end

[ipair,jpair] = find(Cf ~= 0);
MSFpairs = [ipair,jpair]

fprintf(1,'PairIinMSF iEpoch1 iEpoch2\n');
for i=1:numel(ipair)
    fprintf(1,'%5d %5d %5d\n',i,ipair(i),jpair(i));
end

% find simple pair list
PL = pairlist(Qiev);
% print out only those pairs in MSF
fprintf(1,'PairIinMSF Pair iEpoch1 iEpoch2\n');
for i=1:numel(i_msf)
    fprintf(1,'%5d %5d %5d %5d\n',i,i_msf(i),PL(i_msf(i),2),PL(i_msf(i),3));
end
fprintf(1,'The last two columns in the two tables above should match, after sorting.\n');

i_msf
% prune values

TreesMSF = nan(ntrees,nvertices);
if nargout == 6
    QievMSF = Qiev(i_msf,:);
    WeightsMSF = Weights(i_msf);
% Bug here
%     for i=1:ntrees
%         TreesMSF(i,:) = Trees(i,i_msf');
%     end

end


return
end

