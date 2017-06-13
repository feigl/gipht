function [a,b,c,d] = plane_equation(pt1,pt2,pt3)
% given coordinates of three points as column vectors, find equation of
% plane
% https://en.wikipedia.org/wiki/Plane_(geometry)
% Conversely, it is easily shown that if a, b, c and d are constants and a,
% b, and c are not all zero, then the graph of the equation
% 
%         a x + b y + c z + d = 0 , {\displaystyle ax+by+cz+d=0,}
%         ax+by+cz+d=0,
% 
% is a plane having the vector n = (a, b, c) as a normal.[4] This familiar
% equation for a plane is called the general form of the equation of the
% plane.[5]
% 
% Thus for example a regression equation of the form y = d + ax + cz (with
% b = ?1) establishes a best-fit plane in three-dimensional space when
% there are two explanatory variables.


X=[pt1(1),pt2(1),pt3(1)]';
Y=[pt1(2),pt2(2),pt3(2)]';
Z=[pt1(3),pt2(3),pt3(3)]';


%% try again with least squares
c=-1;
G = [X,Y,[1;1;1]]/c;
mest = lscov(G,Z);
a=mest(1);
b=mest(2);
d=mest(3);

 

end

