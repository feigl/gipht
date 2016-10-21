function UHS = wang539(volume_strain,nu,depth,obse,obsn)
%function UHS = wang539(volume_strain,nu,depth,obse,obsn)
% Green's function for a temperature change in a 
% strain nucleus located at (X,Y,Z) = (E, N, Depth) = (0,0,zprime)
% for an observation point located at on the surface at
% (X,Y,Z) = (E, N, Depth) = (obsn, obse, 0)

% find displacement at surface due to volumetric strain
% Wang, H. F. (2000), Theory of poroelasticity with applications to
% geomechanics and hydrology, Princeton University Press.
% 
% Mindlin, R. D., and D. H. Cheng (1950), Nuclei of Strain in the
% Semi-Infinite Solid, Journal of Applied Physics, 21, 926.
% 
% Mindlin, R. D., and D. H. Cheng (1950), Thermoelastic Stress in the
% Semi-Infinite Solid, Journal of Applied Physics, 21, 931.
% 
% HI Kurt ? My Eqn. 5.30, which obsegrates over the
% increment-of-fluid-content distribution function, means that the
% corresponding Green?s function also is the one for
% increment-of-fluid-content. If you divide Mindlin & Cheng?s spherical
% inclusion solution for a thermoelastic half space by the volume of a
% sphere, you get Eqn. 5.30 but with their beta in place of my c_m.  It is
% the Green?s function for a change in temperature at each spatial poobs.
% The spatial part of all of these Green?s functions (in the braces) are
% the same.  To summarize, use Mindlin & Cheng?s beta in place of c_m and
% obsegrate over the temperature distribution instead of the zeta
% distribution. -          Herb
%  
%  
% Thanks for this. I understand that I can use the U^{*HS} vector from
% equation your 5.39 to replace the u_i^* term in the triple obsegral over
% the volume D in your equation 5.30 (page 102 of your book) , but what, if
% anything, should I use for the increment of fluid content
% \zeta(\xi_1,\xi_2,\xi_3)?
%  
%  
% A brief followup. The eqn. for vec{u_e}  in the righthand column of p.
% 932  in the Mindlin & Cheng Thermoelastic paper is the same as Eqn. 5.39
% in my poroelasticity book.  My eqn. is a Green?s function, i.e.,
% displacement for unit volume of temperature change and c_m in my eqn.
% becomes alpha*T*(1+nu)/(1-nu) (Eqn. 1 in Mindlin & Cheng).  So if you
% have a temperature distribution, it should be a straightforward numerical
% obsegration to getting displacement components at the surface. -
% Herb
%  

%% define vertical coordinate of nucleus
zp = depth;

%% define location of nucleus of strain
Xnuc = [0.0; 0.0; zp];

%% define location of image point at depth zp == zprime
Ximg = [0.0; 0.0; -1*zp];

%% define location of observation points 
Xobs = [obsn; obse; 0.];
z = 0.;

%% define vector pointing from nucleus to observation point
R1 = Xobs - Xnuc;

%% define vector pointing from image point to observation point
R2 = Xobs - Ximg;

%% define unit vector pointing downward along positive Z axis
Khat = [0.0; 0.0; 1.0];

%% define beta coefficient by equation (1) of Mindlin and Cheng [1950]
% volume_increase = dT * alpha
beta = volume_strain * (1.0 + nu) / (1.0 - nu);

%% extra factor of 3 to match Okada solution -- Kurt Feigl 20161011
% needed because The factor of 3 in the Okada solution is because the beta
% in Mindlin & Cheng involves linear
% thermal expansion, not volumetric.
beta = beta / 3.;

%% evaluate equation 5.38 (or 5.39?) of Wang [2000], term by term

UHS1 = R1/power(norm(R1),3);

UHS2 = (3.0-4*nu)* R2 / power(norm(R2),3);

UHS3 = -6.0*z*(z+zp) * R2 / power(norm(R2),5);

UHS4 = (-2.0*Khat/power(norm(R2),3)) * ((3.0-4.0*nu)*(z+zp)-z);

UHS = (beta/4.0/pi) * (UHS1 + UHS2 + UHS3 + UHS4);

%% ad-hoc factor of 3 to match Okada solution
%UHS = UHS / 3.;

return





