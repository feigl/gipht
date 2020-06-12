function unit_vector_east = sunitE(track)
% given satellite track id, return east component of unit vector pointing from target on ground to sensor in orbit along line of sight
% suitable for use with COMSOL
% 20191108 Kurt Feigl
switch track
    case 112 % ALOS Track 112
        unit_vector_east = -0.587837519761173;
    otherwise
        warning(sprintf('Unknown track %d. Setting arbitrarily:',track));
        unit_vector_east = sqrt(3)
end
return
end

