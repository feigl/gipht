function [i_msf, i_rep, Nc] = plotMSF( Qiev, trees, rs, scalefactor)
% function [i_msf, i_rep, Nc ] = plotMSF( DD, trees, rs, scalefactor)
% Using the incidence matrix Q, trees, epochs, and rates, returns
% the minimim spanning forest of the full data set using Prim's Algorithm with edges weighted by the
% uncertainty of the data.  In this graph, nodes represent epochs and edges
% represent pairs (interferograms).  
%
% INPUTS:
%   Qiev        - edge-vertex incidence matrix 
%   trees       - matrix of components in data set
%   rs          - uncertainties of pairwise data 
%   scalefactor - scaling value for data set
%
% OUTPUTS:
%   i_msf       - indices of epochs in minimum spanning forest
%   20200326  KLF corrected line above to read as below
%   i_msf       - indices of pairs (edges) in minimum spanning forest - 
%   i_rep       - indices of epochs eliminated from minimum spanning forest
%   20200326  KLF corrected line above to read as below
%   i_rep       - indices of pairs (edges) eliminated from minimum spanning forest
%   Nc          - vector the length of ntrees showing the number of cycles
%                 per each component
%
% Elena C. Baluyut, UW-Madison  2015-02-20
% 20200326 


% Find pairs
pairlst = pairlist(Qiev);
[PN] = pairlist_tree(pairlst, trees, rs/scalefactor);
[ntrees mdummy] = size(trees);
fprintf(1,'PN is matrix containing information of pairs\n');
fprintf(1,'Itree Ipair is it rsofpair\n')
 
PN

% Round weights to nearest 100th
rRSs = rs/scalefactor*100;
rRSs = round(rRSs);
Wrs = rRSs/100;

%Find mimnimum spanning tree (Prim)
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



%% TODO fix mysterious plotting error from  with Matlab 8.5.0.197613 (R2015a) 
% % Error using cell2mat (line 52)
% % CELL2MAT does not support cell arrays containing cell arrays or objects.
% % 
% % Error in biograph.biograph/hgCorrectFontSize>mycell2mat (line 43)
% %     m = cell2mat(c);
% % 
% % Error in biograph.biograph/hgCorrectFontSize (line 34)
% %    set(mycell2mat(get(mycell2mat(get(h.Edges,'hgline')),'UserData')),'FontSize',edgeFontSize)
% 
% 
% Build biograph
if exist('biograph','file') == 2
%if false
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
    set(bg2, 'LayoutType','equilibrium');
    
    % Print biograph to figure with title
    g = biograph.bggui(bg2);
    f = figure();
    copyobj(g.biograph.hgAxes,f);
    v = axis;
    handle=title({'   ';'Minimum spanning forest separated by tree'});
    set(handle,'HorizontalAlignment', 'center');
    set(handle, 'VerticalAlignment', 'top');
else
    warning(sprintf('%s\n%s\n%s\n' ...
        ,'Matlab function named "biograph" is not included in your search path.'...
        ,'This function depends on the BIOINFORMATICS toolbox.'...
        ,'Skipping plot...'));
%     %% TODO plot using functions from graph toolbox
%     % make list of indices
%     s = [1 2 1 3 2 3 3 3]';
%     t = [2 1 3 1 3 4 5 6]';
%     
%     
%     % make graph using matlab functions
%     G = digraph(s,t);
%     
%     plot(G,'Layout','force');
end

return
end

