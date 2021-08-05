function [UTM_easting_in_meters,UTM_northing_in_meters, UTM_height_in_meters] = xyz_porotomo2utm(x_rotated_in_meters, y_rotated_in_meters, z_rotated_in_meters)
%function [UTM_easting_in_meters,UTM_northing_in_meters] = xy_porotomo2utm(x_rotated_in_meters, y_rotated_in_meters)
% given (x,y) in rotated coordinate system, return (e,n) in UTM
% 20170305 Elena Reinisch, based on code from Kurt Feigl, Mike Cardiff, and
% Christina Morency

%% Intialize
%Points used by Morency in convertutm2xy.pl to define angle necessary
%for rotation. The arrow pointing from the westernmost point to the
%southernmost point of the box was used to define the x axis, with the
%westernmost point representing the origin.
i=0;
i=i+1;bradybox_UTM(i,:) = [ 328761.1411         4408840.1323]; % Top Left corner 
i=i+1;bradybox_UTM(i,:) = [ 327850.8122         4407606.2053]; % Top Right
i=i+1;bradybox_UTM(i,:) = [ 328221.6320         4407332.5948];
i=i+1;bradybox_UTM(i,:) = [ 329137.6472         4408559.8842];

% Z-axis
origin_UTM_z = 800; % m
origin_porotomo_z = 0; % m
delta_z = origin_porotomo_z - origin_UTM_z; % m

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

% input data is in PoroTomo x,y coordinates
inlist(:,1) = x_rotated_in_meters;
inlist(:,2) = y_rotated_in_meters;

% rotate x,y coord to UTM using the angle adjacent to theta and taking into
% account origin_UTM
outlist = rotate_2d(inlist, 360 - theta) + repmat(origin_UTM,size(inlist,1),1);

UTM_easting_in_meters = outlist(:,1);
UTM_northing_in_meters = outlist(:,2);
UTM_height_in_meters = z_rotated_in_meters - delta_z;

return

