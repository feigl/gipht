function U =penny(xo,yo,xs,ys,nu,E,H,A,P)
% Calculate displacement vector for penny-shaped crack model

%SUN69   Deformation from penny-shaped crack in elastic half-space.
%  [Ur,Uz] = SUN69(R,H,A,V) or SUN69(R,H,A,P,E,nu) computes radial and
%  vertical displacements Ur and Uz on the free surface, due to a 
%  horizontal circular fracture formed in a semi-infinite elastic medium,
%  with following variables:
%       R: radial distance of observation,
%       H: depth of the center of the source from the surface,
%       A: radius of the source with the hydrostatic pressure,
%       V: volume of injected material,
%       P: change of the hydrostatic pressure in the crack.
%       E: Young's modulus,
%      nu: Poisson's ratio (default is 0.25 for isotropic medium).
%
%  [Ur,Uz,B] = SUN69(...) returns also the maximum separation of fracture
%  (vertical displacement) B.
%
%  Equations from Sun [1969], with approximation H/A >> 1. If H/A > 2, error
%  is about 2 to 3%; if H/A > 5, solution is almost perfect.
%
%  Notes:
%     - Equations are all vectorized, so variables R,H,A,V or P can be 
%       scalar or any of them vector or matrix, then outputs will have 
%       the same size.
%     - Convention: Uz > 0 = UP, f is depth so in -Z direction.
%     - Units should be constistent, e.g.: R, H, A, Ur and Uz in m imply
%       V in m3; optional E and P in Pa, nu dimensionless.

ndata = numel(xo);
if numel(yo) ~= ndata
    error('dimension mismatch');
else
    U = zeros(3,ndata);
end

% explicitly assign source parameters
xn = xo-xs;       % easting  coordinate of observation point w.r.t. source
yn = yo-ys;       % northing coordinate of observation point w.r.t. source
R = sqrt(xn.^2 + yn.^2);  % radial distance
theta = atan2(yn,xn);  % angle from observation point w.r.t. source (radians trigonometric sense CCW from east)

[Ur,Uz] = sun69kf(R,H,A,P,E,nu); % Source is specified in terms of pressure

% Rotate horiz. displacements back to the orig. coordinate system:
U(1,:) =Ur.*cos(theta); % eastward  component of displacement vector
U(2,:) =Ur.*sin(theta); % northward component of displacement vector
U(3,:) =Uz;            % vertical component of displacement vector, positive upward

return




