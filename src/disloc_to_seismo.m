function [easts,norths,depths,names]=disloc_to_seismo(params)
%
%given center of upper edge of fault plane in disloc format
% X disloc
% Y disloc
% L disloc
% W disloc
% dip okada
% strike okada (cw from N)
% depth disloc
% params
%   1       2        3      4       5       6      7       8      9     10
% Length - Width - Depth - Dip - Strike - East - North - Sslip - Dslip - Op
% (km)     (km)    (km)    (deg) (deg)    (km)   (km)    (m)     (m)     (m)
%                  upper   neg-  
%                  edge    ative
%uokada = disloc(params(5:14),xloc,nu);
%
% input is center of upper edge
l = params(1);
w = params(2);
dip = -1*params(4);  % 
%dip = params(4);
strike = params(5);
ec = params(6);
nc = params(7);
hc = params(3);


%Coordinates of centroid wrt lower left
rdip=rad(dip);
xx=l/2;
yy=(w/2)*cos(rdip);
hh=(w/2)*sin(rdip);


%form rotation matrix:             
theta = -1*rad(90-strike); 
r(1,1) = cos(theta);
r(1,2) = sin(theta);
r(2,2) = cos(theta);
r(2,1) = -sin(theta); 

% Find lower left
xv(1) = -xx;  
xv(2) = -2*yy;
%rotation Okada towards UTM
u  = r*xv';
% translate 
e00 = ec - u(1);
n00 = nc - u(2);
h00 = hc + 2*hh;


% Find CENTROID 11
xv(1) = xx;  
xv(2) = yy;
%rotation Okada towards UTM
u  = r*xv';
% translate 
e11 = e00 - u(1);
n11 = n00 - u(2);
h11 = h00 - hh;


% find UPPER RIGHT CORNER 22
xv(1) = 2*xx;  
xv(2) = 2*yy;
%rotation Okada towards UTM
u  = r*xv';
e22 = e00 - u(1);
n22 = n00 - u(2);
h22 = h00 - 2*hh;

% find LOWER RIGHT CORNER 20
xv(1) = 2*xx;
xv(2) = 0.;
u  = r*xv';
e20 = e00 - u(1);
n20 = n00 - u(2);
h20 = h00;

% FIND UPPER LEFT CORNER 02
xv(1) = 0.;
xv(2) = 2*yy;
u  = r*xv';
e02 = e00 - u(1);
n02 = n00 - u(2);
h02 = h00 - 2*hh;

% FIND CENTER OF UPPER EDGE
xv(1) = xx;
xv(2) = 2*yy;
u  = r*xv';
e12 = e00 - u(1);
n12 = n00 - u(2);
h12 = h00 - 2*hh;


easts(1) = e00; norths(1) = n00; depths(1) = h00; names{1} = 'LL';  
easts(2) = e20; norths(2) = n20; depths(2) = h20; names{2} = 'UL';  
easts(3) = e22; norths(3) = n22; depths(3) = h22; names{3} = 'UR';  
easts(4) = e02; norths(4) = n02; depths(4) = h02; names{4} = 'LR';  
easts(5) = e00; norths(5) = n00; depths(5) = h00; names{5} = 'LL';  % repeat LL to close rectangle
easts(6) = e22; norths(6) = n22; depths(6) = h22; names{6} = 'UR';  % diagonal
easts(7) = e12; norths(7) = n12; depths(7) = h12; names{7} = 'UC';  % left half of top edge
easts(8) = e02; norths(8) = n02; depths(7) = h02; names{8} = 'LR';  % right half of top edge
easts(9) = e20; norths(9) = n20; depths(8) = h20; names{9} = 'UL';  % other diagonal
easts(10) = e11;norths(10) = n11;depths(10)= h11; names{10}= 'CC';  % centroid

% for i=1:4
%    fprintf(1,'%12.4f %12.4f %12.4f\n',easts(i),norths(i),depths(i));
% end

return;

