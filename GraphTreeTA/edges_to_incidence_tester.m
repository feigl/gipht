% test routines for finding incidence matrix

clear all
itest = 2;
switch itest
    case 1
        
        % make list of indices
        s = [1 2 1 3 2 3 3 3]';
        t = [2 1 3 1 3 4 5 6]';
        
        
        % find edge-vertex matrix using home-made function
        Qiev = edges_to_incidence([s, t]);
        sortrows(Qiev)
        
        % make graph using matlab functions
        G = digraph(s,t);
        
        % find incidence matrix using Matlab function from Bioinformatics toolbox
        I = transpose(full(incidence(G)));
        sortrows(I)
        figure
        plot(G,'Layout','force');
        
        
        % look for differences
        diff = sortrows(Qiev) - sortrows(I)
        
        find(abs(diff) > eps)
        
        % go back
        st = incidence_to_edges(Qiev)
        
    case 2
        
        % two trees]
        % make list of indices
        s = [1 2 4]';
        t = [2 3 5]';
        iedges = [s, t]
        % find edge-vertex matrix using home-made function
        Qiev = edges_to_incidence(iedges);
        %verbose = 1;
        verbose = 0;
        [itrees,Trees] = find_trees_from_incidence(Qiev,verbose);
        
        ipairs = trees_to_pairs(itrees,iedges)
        
        
end





