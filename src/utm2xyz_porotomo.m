function [x_rotated_in_meters, y_rotated_in_meters, z_rotated_in_meters] = utm2xyz_porotomo(UTM_easting_in_meters,UTM_northing_in_meters, UTM_height_in_meters)
%function [x_rotated_in_meters, y_rotated_in_meters] = utm2xy_porotomo(UTM_easting_in_meters,UTM_northing_in_meters)
% given (e,n) in UTM, return (x,y) in rotated coordinate system
% 20170305 Elena C Reinisch and Kurt Feigl, based on code from Mike Cardiff and Christina
% Morency

%% Intialize
%Points used by Morency in convertutm2xy.pl to define angle necessary
%for rotation. The arrow pointing from the westernmost point to the
%southernmost point of the box was used to define the x axis, with the
%westernmost point representing the origin.
i=0;
i=i+1;bradybox_UTM(i,:) = [ 328761.1411         4408840.1323]; 
i=i+1;bradybox_UTM(i,:) = [ 327850.8122         4407606.2053];
i=i+1;bradybox_UTM(i,:) = [ 328221.6320         4407332.5948];
i=i+1;bradybox_UTM(i,:) = [ 329137.6472         4408559.8842];

% Z-axis
origin_UTM_z = 800; % m
origin_porotomo_z = 0; % m
delta_z = origin_UTM_z - origin_porotomo_z; % m

% X-axis points 
xaxisvec_UTM = bradybox_UTM(3,:) - bradybox_UTM(2,:);

% origin
origin_UTM = bradybox_UTM(2,:);

%Degrees of rotation from Northing/Easting for coordinate rotation.
theta = atan2d(xaxisvec_UTM(2),xaxisvec_UTM(1));

%theta   -36.4219 % degrees

%Based on knowing the origin and rotation angle, create anonymous functions
%for translation to and from UTM coordinates.
% utm2xy = @(inlist,org_UTM,theta) rotate_2d(inlist - repmat(org_UTM,size(inlist,1),1),theta);
% xy2utm = @(inlist,org_UTM,theta) rotate_2d(inlist,-theta) + repmat(org_UTM,size(inlist,1),1);

inlist(:,1) = UTM_easting_in_meters;
inlist(:,2) = UTM_northing_in_meters;

outlist = rotate_2d(inlist - repmat(origin_UTM,size(inlist,1),1),theta);

x_rotated_in_meters = outlist(:,1);
y_rotated_in_meters = outlist(:,2);
z_rotated_in_meters = UTM_height_in_meters - delta_z;

return

