function c=rarcm(a,b)
% given phase gradient values a and b each in radians on [-pi,pi] 
% return arc(a,b) on [0 pi]
% see eq 2.3.13 page 19, Mardia and Jupp [2000]
c = pi - abs(pi-abs(a-b));
return

