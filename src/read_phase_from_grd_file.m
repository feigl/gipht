function [xgrd,ygrd,phaimg] = read_phase_from_grd_file(fn0)
% function [xgrd,ygrd,phaimg] = read_phase_from_grd_file(fn0)
% read phase values from a .grd file written by GMT
% input:
%    fn0 = file name
% output:
%    xgrd = easting in meters or longitude in degrees, unchanged
%    ygrd = northing in meters or latitude in degrees, unchanged
%    phaimg = phase in DN from -128 to +127
% 20130718 Kurt Feigl

fprintf(1,'Reading grid file named %s\n',fn0);
% read grid file
%             D = GRDINFO('file.grd') will in addition return a vector containing
%             (xmin, xmax, ymin, ymax, zmin, zmax, format, xinc, yinc). Format is
%             1 for pixel registration and 0 for grid node registration.
grdinfo = grdinfo2(fn0);

if numel(grdinfo) < 9
    error(sprintf('Not enough information in grid file named %s\n'),fn0);
else
    if grdinfo(7) ~= 1
        error('Found grid node registration. Expected pixel registration.\n');
    else
        [xgrd,ygrd,phaimg]=grdread2(fn0);
        if ceil(grdinfo(6)-grdinfo(5)) == 255
            fprintf(1,'File named %s contains phase values in Digital Numbers between -127 and +127\n',fn0);
        elseif grdinfo(6)-grdinfo(5) <= 2*pi && grdinfo(6)-grdinfo(5) > 1.0
            warning(sprintf('After reading %s, converting radians to DN\n',fn0));
            phaimg = phaimg * 128 / pi ;
        elseif grdinfo(6)-grdinfo(5) < 1.0 && grdinfo(6)-grdinfo(5) >= 0.0
            warning(sprintf('After reading %s, converting cycles to DN\n',fn0));
            phaimg = phaimg * 256.0;
        else
            error(sprintf('File named %s contains values with unknown dimensions\n',fn0));
        end
    end
end

return;


