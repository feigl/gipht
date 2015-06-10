function outfmt = getfmt(value,name,imode)
%function outfmt = getfmt(value,name)
% given a value and a name, return a format statement
% 2013-JAN-15 Kurt Feigl
if exist('imode','var') == 0
    imode = 0;
end


switch imode
    case 0 % long string
        if abs(value) > 0.1 && abs(value) < 10000
            outfmt = '%2s %3d %#28s %10.4f %10.4f %10.4f %10.4f %7.2f %10.4f\n';
            %  elseif numel(strfind(name,'_in_m_')) > 0
        elseif numel(strfind(name,'Easting')) > 0 || numel(strfind(name,'Northing')) > 0
            outfmt = '%2s %3d %#28s %10.1f %10.1f %10.1f %10.1f %7.2f %10.1f\n';
        else
            outfmt = '%2s %3d %#28s %10.2e %10.2e %10.2e %10.2e %7.2g %10.2e\n';
        end
    case 1
        if abs(value) > 0.1 && abs(value) < 10000
            outfmt = '%#28s %10.4f +/- %10.4f\n';
            %  elseif numel(strfind(name,'_in_m_')) > 0
        elseif numel(strfind(name,'Easting')) > 0 || numel(strfind(name,'Northing')) > 0
            outfmt = '%#28s %10.1f +/- %10.1f\n';
        else
            outfmt = '%#28s %10.2e +/- %10.2e\n';
        end
    otherwise
        warning(sprintf('unknown imode = %d\n',imode));
        outfmt = '%f\n';
end

return
