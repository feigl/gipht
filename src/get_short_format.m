function outfmt = get_short_format(value,name)
%function outfmt = get_short_format(value,name)
% given a value and a name, return a format statement
% 2011-OCT-04 Kurt Feigl
if abs(value) > 0.1 && abs(value) < 10000
    outfmt = '%#28s %3d %10.4f +/- %10.4f\n';
    %  elseif numel(strfind(name,'_in_m_')) > 0
elseif numel(strfind(name,'Easting')) > 0 || numel(strfind(name,'Northing')) > 0
    outfmt = '%#28s %3d %10.1f +/- %10.1f\n';
else
    outfmt = '%#28s %3d %10.2e +/- %10.2e\n';
end
return
