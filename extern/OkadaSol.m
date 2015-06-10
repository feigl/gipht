% This file is distributed in the hope that it will be useful, but  
% WITHOUT ANY WARRANTY; without even the implied warranty of  
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% -----------------------------------------------------------------------
% Matlab program for the Okada solution computation
% -----------------------------------------------------------------------
%%% Author: Denys Dutykh, CNRS - University College Dublin
%%% E-mail: Denys.Dutykh@ucd.ie
%%%    URL: http://www.denys-dutykh.com/
% -----------------------------------------------------------------------

function [u v w] = OkadaSol (X, Y, x0, y0, d, U, theta)
% Explanation of arguments :
%   (X, Y) : points where we compute the displacements
% (x0, y0) : left bottom corner of the fault surface
%        d : depth of the fault (at the same corner, normally)
%        u : slip on this fault
%    theta : rake
%    delta : dip angle
%      phi : strike angle

global delta phi L W

% P-wave velocity
vp = 6000;

% S-wave velocity
vs = 3400;

% crust density
rho = 2700;

% we deduce Lame coefficients
mu = vs^2*rho;
lambda = vp^2*rho - 2*mu;

u = zeros (size (X));
v = zeros (size (X));
w = zeros (size (X));

X1 = (X - x0)*cos(phi) + (Y - y0)*sin(phi);
Y1 = -(X - x0)*sin(phi) + (Y - y0)*cos(phi);

P1 = Y1*cos(delta) + d*sin(delta);

u = -U*(cos(theta)*...
    (u1u (X1,P1,Y1,d,delta,lambda,mu) +...
    u1u (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u1u (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u1u (X1-L,P1,Y1,d,delta,lambda,mu)) +...
    sin(theta)*...
    (u2u (X1,P1,Y1,d,delta,lambda,mu) +...
    u2u (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u2u (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u2u (X1-L,P1,Y1,d,delta,lambda,mu)))/(2*pi);

v = -U*(cos(theta)*...
    (u1v (X1,P1,Y1,d,delta,lambda,mu) +...
    u1v (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u1v (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u1v (X1-L,P1,Y1,d,delta,lambda,mu)) +...
    sin(theta)*...
    (u2v (X1,P1,Y1,d,delta,lambda,mu) +...
    u2v (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u2v (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u2v (X1-L,P1,Y1,d,delta,lambda,mu)))/(2*pi);

w = -U*(cos(theta)*...
    (u1w (X1,P1,Y1,d,delta,lambda,mu) +...
    u1w (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u1w (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u1w (X1-L,P1,Y1,d,delta,lambda,mu)) +...
    sin(theta)*...
    (u2w (X1,P1,Y1,d,delta,lambda,mu) +...
    u2w (X1-L,P1-W,Y1,d,delta,lambda,mu) -...
    u2w (X1,P1-W,Y1,d,delta,lambda,mu) -...
    u2w (X1-L,P1,Y1,d,delta,lambda,mu)))/(2*pi);

function res = u1u (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    dbar = eta*sin(delta) - q*cos(delta);
    X = sqrt (xi.^2 + q.^2);
    R = sqrt (xi.^2 + eta.^2 + q.^2);
    I5 = 2*mu/((lambda+mu)*cos(delta))*...
        atan ((eta.*(X+q*cos(delta))+...
        X.*(R+X)*sin(delta))./(xi.*(R+X)*cos(delta)));
    I1 = -mu/(lambda+mu)*xi./((R + dbar)*cos(delta)) -...
        tan(delta)*I5;
    res = xi.*q./(R.*(R+eta)) + atan(xi.*eta./(q.*R)) + I1*sin(delta);

function res = u2u (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    R = sqrt (xi.^2 + eta.^2 + q.^2);
    dbar = eta*sin(delta) - q*cos(delta);
    ybar = eta*cos(delta) + q*sin(delta);
    I4 = mu/((lambda+mu)*cos(delta))*(log(R+dbar) -...
        sin(delta)*log(R+eta));
    I3 = mu/(lambda+mu)*(ybar./(cos(delta)*(R + dbar)) -...
        log(R + eta)) + tan(delta)*I4;
    res = q./R - sin(delta)*cos(delta)*I3;

function res = u1v (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    R = sqrt (xi.^2 + eta.^2 + q.^2);
    dbar = eta*sin(delta) - q*cos(delta);
    ybar = eta*cos(delta) + q*sin(delta);
    I4 = mu/((lambda+mu)*cos(delta))*(log(R+dbar) -...
        sin(delta)*log(R+eta));
    I3 = mu/(lambda+mu)*(ybar./(cos(delta)*(R + dbar)) -...
        log(R + eta)) + tan(delta)*I4;
    I2 = -mu/(lambda+mu)*log(R+eta) - I3;
    res = ybar.*q./(R.*(R+eta)) + cos(delta)*q./(R+eta) + sin(delta)*I2;

function res = u2v (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    dbar = eta*sin(delta) - q*cos(delta);
    ybar = eta*cos(delta) + q*sin(delta);
    X = sqrt (xi.^2 + q.^2);
    R = sqrt (xi.^2 + eta.^2 + q.^2);
    I5 = 2*mu/((lambda+mu)*cos(delta))*...
        atan ((eta.*(X+q*cos(delta))+...
        X.*(R+X)*sin(delta))./(xi.*(R+X)*cos(delta)));
    I1 = -mu/(lambda+mu)*xi./((R + dbar)*cos(delta)) -...
        tan(delta)*I5;
    res = ybar.*q./(R.*(R+xi)) + cos(delta)*atan(xi.*eta./(q.*R)) -...
        sin(delta)*cos(delta)*I1;

function res = u1w (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    eta = eta + eps;
    R = sqrt (xi.^2 + eta.^2 + q.^2);
    dbar = eta*sin(delta) - q*cos(delta);
    I4 = mu/((lambda+mu)*cos(delta))*(log(R+dbar) -...
        sin(delta)*log(R+eta));
    res = dbar.*q./(R.*(R+eta)) + q*sin(delta)./(R+eta) +...
        I4*sin(delta);
    
function res = u2w (xi, eta, y, d, delta, lambda, mu)
    q = y*sin(delta) - d*cos(delta);
    X = sqrt (xi.^2 + q.^2);
    R = sqrt (X.^2 + eta.^2);
    dbar = eta*sin(delta) - q*cos(delta);
    I5 = 2*mu/((lambda+mu)*cos(delta))*...
        atan ((eta.*(X+q*cos(delta))+...
        X.*(R+X)*sin(delta))./(xi.*(R+X)*cos(delta)));
    res = dbar.*q./(R.*(R+xi)) +...
        sin(delta)*atan(xi.*eta./(q.*R)) -...
        I5*sin(delta)*cos(delta);
