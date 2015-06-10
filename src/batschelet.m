function A=batschelet(kappa)
% given kappa concentration parameter in Von Mises distribution,
% return the mean resultant length A = Rbar
% from Mardia 1972 equation (5.4.5) page 122
% also Mardia and Jupp equation (3.5.31) page 40
% also N.I. Fisher [1993] page 49
A = besseli(1,kappa)/besseli(0,kappa);
return




