function theta = incidence_angle(unitv_east,unitv_north,unitv_up)
%function theta = incidence_angle(unitv_east,unitv_north,unitv_up)
% given unit vector [E,N,U], return incidence angle in degrees from vertical
%
% Example 1: from ERS, using known unit vector
%     incidence_angle(0.3973,-0.1088,0.9112)
%     ans =
%        24.3264
%    
% Example 2: from ERS, using orbit file to calculate unit vector
%     incidence_angle(unitv_east,unitv_north,unitv_up)
%     ans =
%        24.3267
%     [unitv_east,unitv_north,unitv_up] = lookvector(-17.97459677,64.09867826,200,'11176.orb')
%     Read    80 records from orbit file 11176.orb
%     Returning    80 valid records from orbit file 11176.orb
%     unitv_east =
%         0.3973
%     unitv_north =
%        -0.1088
%     unitv_up =
%         0.9112
%     incidence_angle(unitv_east,unitv_north,unitv_up)
%     ans =
%        24.3267
%
%
% 20130629 Kurt Feigl

theta = rad2deg(atan2(hypot(unitv_east,unitv_north),unitv_up));
return

