%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the data structure needed for an ellipsoid inhomogeniety problem;
% solve the problem with an equevilent Eshelby inclusion problem;
% demonstrate the outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Sharma's ellpdoidal heterogeneity with observation line
%clear all;
%poisson ratio of matrix
incl.vm= 0.25;
%elastic modulus of matrix
incl.Em=2.2e10;
%hetergeneity poisson ratio
incl.vh=0;
%hetergeneity elastic modulus
incl.Eh=0;

% dimensiona of the ellipsoid.
incl.dim=[1 .2 20];
% ortation angles in radian applied to the ellipsoid in sequence of [x y z];.
incl.ang = [0 0 0];

% remote stress ordered by sigma11, sigma12, sigma13, sigma22,sigma23,sigma33
incl.stressvec=[-1e6;0;0;0;1e6;0];
incl.eigp=[0;0;0;0;0;0];
% first observaton grid
x = 0;
y = -25:.5:25;
z = -25:.5:25;
incl.grid{1} = {x,y,z};
% second observation grid
x = -25:.5:25;
y = 0;
z = -25:.5:25;
incl.grid{2} = {x,y,z};
% third observation grid
x = -25:.5:25;
y = -25:.5:25;
z = 0;
incl.grid{3} = {x,y,z};
clear x y z
% add more grids if needed

% call Eshelby solver,arg: 'disp','stress','strain' only output displacements by default
incl.sol = Esh_sol(incl,'disp','stress','strain');

% demonstrate the third grid's results for the first and second components of
% the output displacement and stress.
X = squeeze(incl.sol.grid{1,1});
Y = squeeze(incl.sol.grid{1,2});
Z = squeeze(incl.sol.grid{1,3});
u = squeeze(incl.sol.u{1});
stress = squeeze(incl.sol.stress{1});
figure
az = 0;
el = 90;
subplot(2,2,1)
surf(X,Y,stress(:,:,1),'LineStyle','none');
colorbar;
view(az, el);
subplot(2,2,2)
surf(X,Y,stress(:,:,2),'LineStyle','none');
colorbar;
view(az, el);
subplot(2,2,3)
surf(X,Y,stress(:,:,3),'LineStyle','none');
colorbar;
view(az, el);
subplot(2,2,4)
surf(X,Y,stress(:,:,4),'LineStyle','none');
view(az, el);
colorbar;