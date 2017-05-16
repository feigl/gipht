%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ctensord.m
% this function calculates the 4th order elastic stiffness tensor given elastic
% modulus and poisson's ratio for an isotropic material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function q=Ctensord(Em,vm)


Gm=Em/(2+2*vm);
lamem=2*Gm*vm/(1-2*vm);
q=zeros(6,6);

q=[lamem+2*Gm,0,0,lamem,0,lamem;
0,2*Gm,0,0,0,0;
0,0,2*Gm,0,0,0;
lamem,0,0,lamem+2*Gm,0,lamem;
0,0,0,0,2*Gm,0;
lamem,0,0,lamem,0,lamem+2*Gm];