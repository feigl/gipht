function [i_msf, i_rep, Nc ] = plotMSF( DD, trees, rs, scalefactor)
% function [i_msf, i_rep, Nc ] = plotMSF( DD, trees, rs, scalefactor)
% Using the incidence matrix Q, trees, epochs, and rates, returns
% the minimim spanning forest of the full data set using Prim's Algorithm with edges weighted by the
% uncertainty of the data.  In this graph, nodes represent epochs and edges
% represent pairs (interferograms).  
%
% INPUTS:
%   Q           - edge-vertex incidence matrix 
%   trees       - matrix of components in data set
%   rs          - uncertainties of pairwise data 
%   scalefactor - scaling value for data set
%
% OUTPUTS:
%   i_msf       - indices of epochs in minimum spanning forest
%   i_rep       - indices of epochs eliminated from minimum spanning forest
%   Nc          - vector the length of ntrees showing the number of cycles
%                 per each component
%
% Elena C. Baluyut, UW-Madison
% 2015-02-20

% Find pairs
pairlst = pairlist(DD);
[PN] = pairlist_tree(pairlst, trees, rs/scalefactor);
[ntrees mdummy] = size(trees);

% Round weights to nearest 100th
rRSs = rs/scalefactor*100;
rRSs = round(rRSs);
Wrs = rRSs/100;

%Find and plot mimnimum spanning tree (Prim)
[ i_msf, i_rep, Nc ] = MSF(ntrees, trees, PN);
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
    Cf(PN_msf(h,3), PN_msf(h,4)) = rs(h)/scalefactor; % store uncertainty for each pair in MSF
end

% Build biograph
bgf = biograph(Cf);
ids = get(bgf.nodes,'ID');

% Add labels and assign weights
for h = 1:mn
    node = num2str(h);
    ids{h} = ['epoch', ' ', node]; % add epoch labels to biograph
end
bg2 = biograph(Cf, ids);
set(bg2, 'ShowArrows', 'on') % DAG representation (shows chronological order of epochs in pair)
set(bg2, 'ShowWeights', 'off') % don't show edge weights

% Print biograph to figure with title
g = biograph.bggui(bg2);
f = figure();
copyobj(g.biograph.hgAxes,f);
v = axis
handle=title({'   ';'Minimum spanning forest separated by tree'});
set(handle,'HorizontalAlignment', 'center');
set(handle, 'VerticalAlignment', 'top');
end

