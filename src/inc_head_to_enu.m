function [ue, un, uz] = inc_head_to_enu(incidence,heading)
%function [ue, un, uz] = inc_head_to_enu(incidence,heading)
% calculate unit vector s^ pointing from ground to satellite
% inputs:
%   incidence == incidence angle in degrees from vertical
%    heading   == azimuth of ground track of sensor in degrees clockwise from north
% outputs:
%     ue == east   component of unit vector
%     un == north  component of unit vector
%     uz == upward component of unit vector
% assumptions:
%     sensor looks right
%     
%     20150615 Kurt Feigl
     
%%     an example by Helene Le Mevel
% hengill.geology.wisc.edu% more /data/chile/UAVSAR/calc_unitvector.m
% % UAVSAR
% %pair 1 DESC
% % cos is in radian in Matlab!!
% inc=53.64*pi/180;
% %az=(286.412-180);
% az=286.412;
% theta=(az-90)*pi/180;
% 
% 
%%CSK DESCENDING
% inc=45.2237548*pi/180;
% az=(90-77.077972412)+180;
% theta=(az-90)*pi/180;
% % 
% ue=sin(inc)*sin(theta)
% un=sin(inc)*cos(theta)
% uz=cos(inc)
% 
% [unitv_e,unitv_n,unitv_u] = inc_head_to_enu(45.22375,(90-77.077972412)+180)
% ue =
%     0.6919
% un =
%    -0.1587
% uz =
%     0.7043
% unitv_e =
%     0.6919
% unitv_n =
%    -0.1587
% unitv_u =
%     0.7043

% 
%% CSK pair2 ASC
% inc2=44.9884*pi/180;
% az2=106.339;
% theta2=(az2-90)*pi/180;
% 
% ue2=sin(inc2)*sin(theta2)
% un2=sin(inc2)*cos(theta2)
% uz2=cos(inc2)

% [unitv_e,unitv_n,unitv_u] = inc_head_to_enu(44.98,106)
% ue2 =
%     0.1989
% un2 =
%     0.6784
% uz2 =
%     0.7072
% unitv_az =
%     16
% unitv_e =
%     0.1948
% unitv_n =
%     0.6795
% unitv_u =
%     0.7074


% sensor looks to the right, i.e. more clockwise
look_az = heading + 90.; 

% unit vector points in opposite direction
unitv_az = mod(look_az + 180., 360.);

% trigonometric angle measured COUNTER clockwise from X (easting) axis
theta = (90. - unitv_az) * pi / 180.;

inc = incidence * pi/180;

ue=sin(inc)*cos(theta);
un=sin(inc)*sin(theta);
uz=cos(inc);


end

