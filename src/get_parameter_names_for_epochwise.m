function [pnames,pscl] = get_parameter_names_for_epochwise(me)
% function [pnames,pscl] = get_parameter_names_for_epochwise(me)
% return the index and name of parameter pname1
% also return scale factor 'pscl' in same units as parameter.
% The scale factor will be used for numerical calculation of partial
% derivatives. It should be of the same order of magnitude as a
% "reasonable-sized" adjustment to the parameter. 
% For example, the scale factor for a UTM coordinate might be 1000 m.
 
j=0;
for i=1:me
    j=j+1;
    pnames{i} = sprintf('time fn @ epoch %03d in years    ',j);
    pscl(i)   = 1.0; % years
end
j=0;
for i=1*me+1:2*me
    j=j+1;
    pnames{i} = sprintf('E grad  @ epoch %03d dimless     ',j);
    pscl(i)   = 1.0e-6; % 1 mm over 1 km    
end
j=0;
for i=2*me+1:3*me
    j=j+1;
    pnames{i} = sprintf('N grad  @ epoch %03d dimless     ',j);
    pscl(i)   = 1.0e-6; % 1 mm over 1 km    
end
j=0;
for i=3*me+1:4*me
    j=j+1;
    pnames{i} = sprintf('U grad  @ epoch %03d dimless     ',j);
    pscl(i)   = 1.0e-6; % 1 mm over 1 km    
end
j=0;
for i=4*me+1:5*me
    j=j+1;
    pnames{i} = sprintf('OrbitHoriz @ epoch %03d in m     ',j);
    pscl(i)   = 1.0; % meters
end
j=0;
for i=5*me+1:6*me
    j=j+1;
    pnames{i} = sprintf('OrbitAlong @ epoch %03d in m     ',j);
    pscl(i)   = 1.0; % meters
end
j=0;
for i=6*me+1:7*me
    j=j+1;
    pnames{i} = sprintf('OrbitVerti @ epoch %03d in m     ',j);
    pscl(i)   = 1.0; % meters
end
j=0;
for i=7*me+1:8*me
    j=j+1;
    pnames{i} = sprintf('OrbitVelH  @ epoch %03d m per s  ',j);
    pscl(i)   = 1.0; % m/s
end
j=0;
for i=8*me+1:9*me
    j=j+1;
    pnames{i} = sprintf('OrbitVelA  @ epoch %03d m per s  ',j);
    pscl(i)   = 1.0; % m/s
end
j=0;
for i=9*me+1:10*me
    j=j+1;
    pnames{i} = sprintf('OrbitVelV  @ epoch %03d m per s  ',j);
    pscl(i)   = 1.0; % meters
end
j=0;
for i=10*me+1:11*me
    j=j+1;
    pnames{i} = sprintf('Offset  @ epoch %03d in cycles   ',j);
    pscl(i)   = 0.1; % cycles
end
%     j=0;
%     for i=6*me+1:7*me
%         j=j+1;
%         pnames{i} = sprintf('PixDX   @ epoch %03d in meters   ',j);
%     end
%     j=0;
%     for i=7*me+1:8*me
%         j=j+1;
%         pnames{i} = sprintf('PixDY   @ epoch %03d in meters   ',j);
%     end

pnames = truncate_parameter_names(pnames);
return






