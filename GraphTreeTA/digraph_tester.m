%% make a graph using Matlab functions
clear all
close all
s = [1 1 1 2 2 3 3 4 5 5 6 7];
t = [2 4 8 3 7 4 6 5 6 8 7 8];
weights = [10 10 1 10 1 10 1 1 12 12 12 12];
names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
G = digraph(s,t,weights,names)
figure;
plot(G,'Layout','force','EdgeLabel',G.Edges.Weight)


