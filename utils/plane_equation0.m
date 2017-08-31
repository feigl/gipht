function [a,b,c,pt0] = plane_equation(pt1,pt2,pt3)
% given coordinates of three points as column vectors, find equation of
% plane

X=[pt1(1),pt2(1),pt3(1)]';
Y=[pt1(2),pt2(2),pt3(2)]';
Z=[pt1(3),pt2(3),pt3(3)]';

%% try with least squares
% Xm = mean(X);
% Ym = mean(Y);
% Zm = mean(Z);
% 
% G = [X-Xm,Y-Ym,Z-Zm]
% mest = lscov(G,zeros(3,1));
% pt0 = [Xm,Ym,Zm]';
% a = mest(1);
% b = mest(2);
% c = mest(3);


%% try again with least squares
G = [X,Y,[-1;-1;-1]];
mest = lscov(G,Z);
a=mest(1);
b=mest(2);
c=

% %% try as eigenvalue problem
% Xm = mean(X);
% Ym = mean(Y);
% Zm = mean(Z);
% 
% R = [X-Xm, Y-Ym, Z-Zm];
% %Computation of the principal directions
% %[V,D] = eig(R'*R);
% [V,D] = eig(R);
% 
% %Extract the output from the eigenvectors
% n = V(:,1);
% mest = n;
% V = V(:,2:end);
% pt0 = [Xm,Ym,Zm]';


%% Try with cross product
% v12=pt2-pt1;%v12=v12/norm(v12)
% v13=pt3-pt1;%v13=v13/norm(v13)
% mest=cross(v12,v13)/norm(v12)/norm(v13);
%mest=mest/norm(mest);
%pt0 = pt1;

 

end

