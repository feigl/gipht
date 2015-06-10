function best_tour = tspsiman(EUC_2D,Tend) 
% TSPIMAN A Symmetric 2D Euclidian Traveling Salesman Problem (TSP) 
% Nearest Neighbor tour construction + 2-Opt local search + Simulated Annealing/Metropolis test
% Developed using MATLAB v5.
%
% Main program by Jericho A. Corpus - Dept. of Mathematics, Graduate Studies 
%  University of the Philippines, Diliman, Manila - March 1998
%  
% For inquiries, comments, corrections, e-mail me: drcorps@hotmail.com
%
% slightly modified by Kurt Feigl 2004 JUL 11
%
% This function makes use of the following custom function/s and data file/s:
%  tourdist.m 
%  kroa100.m - 100 cities in Euclidian 2D format, best recorded lowest bound = 21282
%
% To run program, type the m-filename of the modified TSP data file in the MATLAB prompt
%  (in this case type "kroa100") and ENTER.
%
% A. Input
% Sample instances from the Travelling Salesman Library (TSPLIB) in
%  ftp://softlib.es.rice.edu/pub/tsplib, or
%  http://www.iwr.uni-heidelberg.de/iwr/comopt/soft/TSPLIB95/TSPLIB.html
%
% Weights or lengths must be given in x-y coordinates (EUC_2D). Data file must be slightly 
%  modified and saved as a regular text file for it to be processed by the MATLAB program. 
%  See kroa100.m (KroA100 in TSPLIB).
%
% B. Output
%  Best_tour_length
%  Temperature_of_best_tour_length
%  Solution_count
%  Search_stop_temperature  
%  Elapsed_time % In seconds
%  Solutions_generated - number of solutions generated
%  Floating_point_operations
%  Graphs
%   Temperature vs. Tour Length
%   Number of solutions vs. Tour Length
%
% Some References
% 
% Charon, I. and O. Hudry (1995). Mixing Different Components of Metaheuristics, 
%  Heuristics, 35, 589-603.
%
% Kirkpatrick, S., C.D. Gelatt, and P.M. Vecchi (1983). Optimization by Simulated 
%  Annealing, Science, 220, 670-680.
%
% Metropolis, W., A. Rosenbluth, M. Rosenbluth, A. Teller, and E.Teller (1953). Equation
%  of the State Calculations by Fast Computing Machines, Journal of Chemical Physics, 21,
%  1087-1208.
%
% Yaguira, M. and T. Ibaraki (1995). Genetic and Local Search Algorithms as Robust and
%  Simple Optimization Tools, Heuristics, 5, 63-82.
besttour = NaN;  % KF initialize
%flops(0) ;
t0= clock ;
xy= EUC_2D ; 
n_cities= length(xy) ;
rand('state',sum(100*clock));

% Create the distance matrix for all of the cities given x-y coordinates
distance_matrix = zeros(n_cities) ;
for n_cities_x = 1: n_cities,
     for n_cities_y = 1:n_cities_x
          x = xy(n_cities_x, 1) ;
          y = xy(n_cities_x, 2) ;
          xx = xy(n_cities_y, 1) ;
          yy = xy(n_cities_y, 2) ;
          distance_matrix(n_cities_x, n_cities_y)= ceil(sqrt((x - xx)^2 + (y - yy)^2)) ;
          distance_matrix(n_cities_y, n_cities_x)= distance_matrix(n_cities_x, n_cities_y) ;
    end
end 
% End of matrix construction

% Construct an initial tour using Nearest Neighbor heuristic 
lenbestNN= inf ;
pbestNN= [] ;

% A random selection of the starting city
%prand= randperm(n_cities) ; 
%f=find(prand==1) ;
%prand(f)= prand(1) ;
%p= [1 prand(2)] ;
%i= prand(3) ;

% start at the first city
p = [1 2];
i=1;

% start at the smallest y value
%j = find(xy == min(xy(:,2)))-1;
%p = [1 j];
%i = j-1;
                       
count= 3 ;
    
while count <= n_cities
        NNdist= inf ;
        pp= i ;
        for j= 1: n_cities
             if (distance_matrix(i, j) < NNdist) & (j~=i) & ((j~=p) == ones(1,length(p)))
                NNdist= distance_matrix(i, j) ; 
                pp= j ;
             end           
        end
        p= [p pp ] ; 
        i= pp ;
        count= count + 1 ;
end

% Computing tour cost or length using Tourdist.m function
len= tourdist(p, distance_matrix) ;
 
if len < lenbestNN
      lenbestNN= len ; 
      pbestNN= p ; 
end 
% End of initial tour construction

solnn= [] ;
lenn= [] ; temp= [] ;
soln= 1 ;
% ========================
% A 2-Opt local search
% ========================
lencurr= lenbestNN; 
Best_tour_length= lenbestNN; 
pcurr= pbestNN ; 
pbest= pbestNN ; 
 
% ========================
% Temperature control
% ========================
restart= 0 ; % KURT CHANGED FROM 0
Tstart= 30 ; % Start temperature
%Tend= 1 ; % Stop temperature
%Tend= 4 ; % Stop temperature
Tred= 0.97 ;
T= Tstart ;
Nochange= 2 ; % If after Nochange neighborhood searches, no improvements or 
              %  changes in tour search, annealing complete, break search.
% ========================

lenn= [lenn lencurr] ;
temp= [temp T] ;
solnn= [solnn soln] ;

bb= 0 ;  

while T >= Tend  
   big= n_cities - 1 ; 
   while big >= 3   
      small= big - 2 ;
      while small >= 1
       
         curropt= distance_matrix(pcurr(big),pcurr(big+1)) + distance_matrix(pcurr(small),pcurr(small+1)) ;  
          swap2= distance_matrix(pcurr(small),pcurr(big)) + distance_matrix(pcurr(small+1),pcurr(big+1)) ; 
  
          soln= soln + 1 ;
       
         if swap2 < curropt
            order2= 1: n_cities ;  
            order2=[1:small big:-1:small+1 big+1:n_cities] ; 
 
            pcurr= pcurr(order2) ;
            lencurr= tourdist(pcurr, distance_matrix) ;
             
            lenn= [lenn lencurr] ; 
            temp= [temp T] ; 
            solnn= [solnn soln] ; 
 
            if lencurr < Best_tour_length
	       Best_tour_length= lencurr;       
               pbest= pcurr ;  
               Temperature_of_best_tour_length= T;
               %fprintf (1,'Temperature of best tour_length = %20.2f\r',T); 
               Solution_count= soln;  
               T= Tred * T ; 
               if T <= 3
                  T= 50 ;
               end
            end
            Tcurr= T ;      
            bb= 0 ;
            big= n_cities - 1 ; 
            small= big - 1 ;
            if T <= 3
               T= 10 ;
            end
 
            if T <= Tend ;
               big= 2.9 ; 
               break
             end   
      
         elseif swap2 > curropt
            %r= abs(randn) ;
            r= rand; % where r ranges from 0.0 to 1.0 
            diff= swap2 - curropt ;
            %if r < exp(-(diff) / T)
              if r <= exp(-(diff) / T)
               order2= 1: n_cities ;
               order2=[1:small big:-1:small+1 big+1:n_cities] ; 
                pcurr= pcurr(order2) ;
                    lencurr= tourdist(pcurr, distance_matrix) ;
               
               T= Tred * T ; 
                 bb= 0 ;  
    
            end
          end
           small= small - 1 ;
      end
      big= big - 1 ;   
   end
   bb= bb + 1 ;   
   if T <= Tend | bb > Nochange 
    
%     clc
      Best_tour_length;
   besttour= [pbest -pbest(1)];  
%     Temperature_of_best_tour_length
%      Solution_count
      Search_stop_temperature= T;
      Elapsed_time= etime(clock, t0); % In seconds
      Solutions_generated= soln ;
%     Floating_point_operations = flops
      if bb > Nochange
         No_change= bb;
      end
    end
end      
% End of local search
 
best_tour = besttour;

