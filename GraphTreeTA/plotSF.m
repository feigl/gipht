function [ output_args ] = plotSF( DD, trees, tu, rs, scalefactor )
%plotSF Using the incidence matrix DD, trees, epochs, and rates, returns
%the connectivity graph of the full data set with edges weighted by the
%uncertainty of the data.  In this graph, nodes represent epochs and edges
%represent pairs (interferograms).  This function is used in conjunction
%with plotMSF to find the optimal abridged data set according to Prim's
%algorithm.

% Find pairs
pairlst = pairlist(DD);
[PN] = pairlist_tree(pairlst, trees, rs/scalefactor);
[ndat mdummy] = size(DD);

%round weights to nearest 100th
rRSs = rs/scalefactor*100;
rRSs = round(rRSs);
Wrs = rRSs/100;

% Build foundation of connectivity matrix
C = zeros(numel(tu), numel(tu));
mn = numel(tu);
for h = 1:ndat
    C(PN(h,3), PN(h,4)) = Wrs(h);
end

% Plot connectivity graph
bg = biograph(C);
ids = get(bg.nodes,'ID');

% Add weights
for h = 1:numel(tu)
    node = num2str(h);
    ids{h} = ['epoch', ' ', node];
end

bg1 = biograph(C, ids);
set(bg1, 'ShowArrows', 'off')
set(bg1, 'ShowWeights', 'on')

% Print biograph to figure with title
g = biograph.bggui(bg1);
f = figure();
copyobj(g.biograph.hgAxes,f);
v = axis
handle=title({'   ';'Connectivity graph of full data set'});
set(handle,'HorizontalAlignment', 'center');
set(handle, 'VerticalAlignment', 'top');
end

