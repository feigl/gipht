function [rot_points] = rotate_2d(points,cwd_angle)

%rotate_2d: Code for rotating a set of points in a 2D plane
%
% © 2005-2013 Michael Cardiff, All Rights Reserved.
%
%[rot_points] = rotate_2d(points,cwd_angle)
%where:
%   -rot_points (numpoints x 2) are the coordinates of the rotated points
%   -points (numpoints x 2) are the coordinates of the unrotated points
%   -cwd_angle is the clockwise angle, in degrees, about which the points
%   are rotated



rot_mat = [cosd(-cwd_angle) -sind(-cwd_angle); ...
    sind(-cwd_angle) cosd(-cwd_angle)];
    
rot_points = (rot_mat*(points'))';

