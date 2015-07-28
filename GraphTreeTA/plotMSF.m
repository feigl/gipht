function [i_msf, i_rep, Nc ] = plotMSF( DD, trees, rs, scalefactor)
%plotMSF Using the incidence matrix DD, trees, epochs, and rates, returns
%the minimim spanning forest of the full data set using Prim's Algorithm with edges weighted by the
%uncertainty of the data.  In this graph, nodes represent epochs and edges
%represent pairs (interferograms).  

% Find pairs
pairlst = pairlist(DD);
[PN] = pairlist_tree(pairlst, trees, rs/scalefactor);
[ntrees mdummy] = size(trees);

%round weights to nearest 100th
rRSs = rs/scalefactor*100;
rRSs = round(rRSs);
Wrs = rRSs/100;

%Find and plot mimnimum spanning tree (Prim)
[ i_msf, i_rep, Nc ] = MSF(ntrees, trees, PN);
i_r = [];
for jj = 1:numel(i_rep)
    [i_n mdum] = find(PN(:,2) == i_rep(jj));
    i_r(end+1,:) = i_n;
end
PN_msf = PN;
PN_msf(i_r, :) = [];
max_end = max([max(PN_msf(:,3)) max(PN_msf(:,4))]);

Cf = zeros(max_end, max_end);
mn = max_end;
[hh dum] = size(PN_msf);
for h = 1:hh
    Cf(PN_msf(h,3), PN_msf(h,4)) = rs(h)/scalefactor;
end

bgf = biograph(Cf);
ids = get(bgf.nodes,'ID');

% Add labels and assign weights
for h = 1:mn
    node = num2str(h);
    ids{h} = ['epoch', ' ', node];
end
bg2 = biograph(Cf, ids);
set(bg2, 'ShowArrows', 'on')
set(bg2, 'ShowWeights', 'off')

% Print biograph to figure with title
g = biograph.bggui(bg2);
f = figure();
copyobj(g.biograph.hgAxes,f);
v = axis
handle=title({'   ';'Directed Acyclic Graph Representing Optimal Abridged Data Set'});
set(handle,'HorizontalAlignment', 'center');
set(handle, 'VerticalAlignment', 'top');
end

