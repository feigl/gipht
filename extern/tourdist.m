function tour_distance = tourdist(tourvec, distance_matrix);
% TOURDIST: This function is used with Tspsiman.m.
% 2D Euclidian Traveling Salesman Problem (TSP)
% Nearest Neighbor tour construction + 2-Opt local search + Simulated Annealing/Metropolis test
% Developed using MATLAB v5.
%
% Main TSP program by Jericho A. Corpus - Dept. of Mathematics, Graduate Studies 
%  University of the Philippines, Diliman, Manila - March 1998
%   
% For inquiries, comments, corrections, e-mail me: drcorps@hotmail.com

n_cities = length(tourvec);
city = 1;
tour_distance = 0;

while city <= (n_cities - 1),
   tour_distance = tour_distance + distance_matrix(tourvec(city), tourvec(city+1));
   city = city + 1;
end;

tour_distance = tour_distance + distance_matrix(tourvec(n_cities), tourvec(1));

