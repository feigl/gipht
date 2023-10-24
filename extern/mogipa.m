function u=mogipa(volgeom, xloc, nu, mu)
%function u=mogipa(volgeom, xloc, nu, mu)
%
%Computes surface displacements, strains, and tilts due to a Mogi source
%
%
%Inputs:
%     volgeom = Mogi source geometry: East, North, Depth, Radius, pressure
%               <length, length, length, radius, Pressure>
%        xloc = Matrix of local station coordinates <length>, stored columnwise,
%               i.e., east in first row, north in second row
%          nu = Poisson's ratio
%          mu = shear modulus 
%
%Output:
%           u = matrix of displacements: Ux, Uy, Uz <length>
%           e = matrix of strains: Exx, Exy, Eyy
%           t = matrix of tilts: dUz/dx, dUz/dy
%
%Notes: The term 'depth' denotes an unsigned length and should therefore always be
%given positive.  Keep your length units consistent! If you mix km and m you may get
%unexpected results, particularly with the strains and tilts.   
%
%05-17-98 Peter Cervelli.
% 2023/02/16 Kurt Feigl
% parameterized in terms of Pressure and chamber radius  as 

%Check arguments

	if nargin < 3 | nargin > 5 | nargout > 3
		error('Usage incorrect.');
	end

	[v]=size(volgeom);

	if min(v)~=1 | max(v)~=5
		error('First argument must be a 5 element source geometry vector.');
	end

	[x]=size(xloc);
	
	if x(1)~=2
		error('Second argument must be a 2xn matrix of station coordinates.');
	end

    if length(nu) ~= 1
        error('Fourth argument must be a scalar (dimensionless Poisson''s ratio');
    end

    if length(mu) ~= 1
        error('Fifth argument must be a scalar (shear modulus in Pa).');
    end

%Compute displacements
% Use equations 7.12 from 
% Segall, P. (2010), Earthquake and volcano deformation, 432 pp., Princeton University Press, Princeton, N.J.  


	E=volgeom(1)-xloc(1,:);
	N=volgeom(2)-xloc(2,:);
	E2=E.^2;
	N2=N.^2;
	d2=volgeom(3)^2;
    C = (1-nu)*volgeom(5)*(volgeom(4)^3)/mu;
	R2=(d2+E2+N2);
	R3=C ./ (R2.^(3/2));
	u=[E.*R3; N.*R3; volgeom(3)*R3];

% %Compute strains (if necessary)
% 
% 	if nargout > 1
%  		R5=C*R.^-5;
%      	e(1,:)=R5.*(2*E2-N2-d2);
%      	e(2,:)=3*R5.*(E.*N);
%      	e(3,:)=R5.*(2*N2-E2-d2);
% 	end
% 
% %Compute tilts (if necessary)
% 
% 	if nargout > 2
%        	t(1,:)=3*R5.*E;
%        	t(2,:)=3*R5.*N;
% 	end
