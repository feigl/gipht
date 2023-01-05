function INFO = grdinfo4(file)
%GRDINFO2  Print information about a GMT grdfile (netCDF format, GMT v3 or v4)
%
% Uses NetCDF libraries to display information about a GMT grid file.
% Duplicates (some) functionality of the program grdinfo (which requires
% compilation as a mexcdf function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% GRDINFO2('file.grd') will display information about the GMT grid
% file 'file.grd' in a format similar to the gmt command grdinfo.
%
% D = GRDINFO('file.grd') will in addition return a structure containing
% (xmin, xmax, ymin, ymax, zmin, zmax, format, xinc, yinc). Format is
% 1 for pixel registration and 0 for grid node registration.
%
% See also GRDREAD3, GRDWRITE3

% This program is expected to work on any GMT netCDF format file,
% but it does not duplicate all the functionality of GMT's I/O
% library, so it will not work on all files supported by GMT.
% In particular, it will fail on binary format grdfiles.
% It is the responsibility of the user to determine whether this
% program is appropriate for any given task.
%
% For more information on GMT grid file formats, see:
% http://www.soest.hawaii.edu/gmt/gmt/doc/gmt/html/GMT_Docs/node70.html
% Details on Matlab's native netCDF capabilities are at:
% http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html
%
% GMT (Generic Mapping Tools, <http://gmt.soest.hawaii.edu>)
% was developed by Paul Wessel and Walter H. F. Smith
%
% Kelsey Jordahl
% Marymount Manhattan College
% http://marymount.mmm.edu/faculty/kjordahl/software.html
%
% Time-stamp: <Wed Jan  6 16:26:46 EST 2010>
%
% Version 1.1.1, 6-Jan-2010
% first released on MATLAB Central
% revised 20160813 Kurt Feigl
% revised 20220823

if nargin < 1
    help(mfilename);
    return
end

% check for appropriate Matlab version (>=7.7)
V=regexp(version,'[ \.]','split');
if (str2num(V{1})<7) || (str2num(V{1})==7 && str2num(V{2})<7)
    ver
    error('grdread3: Requires Matlab R2008b or later!');
elseif (str2num(V{1}) >= 9)
    ver9=true;
else
    ver9=false;
end

if ver9 == true
    II=ncinfo(file)

    INFO.source = II.Attributes(2).Value
    INFO.title=file;
    INFO.description='description';

    fprintf(1,'Variables\n');
    II.Variables.Name
    nvars = numel(II.Variables)
%     name1=II.Variables(1).Attributes(1).Value
%     name2=II.Variables(2).Attributes(1).Value
%     name3=II.Variables(3).Attributes(1).Value
%     name4=II.Variables(4).Attributes(1).Value

    if II.Variables(2).Attributes.Name == 'units'
        INFO.xname=II.Variables(2).Attributes.Value
    else
        INFO.xname=' ';
    end
    if II.Variables(3).Attributes.Name == 'units'
        INFO.yname=II.Variables(3).Attributes.Value
    else
        INFO.yname=' ';
    end
    if II.Variables(4).Attributes.Name == 'units'
        INFO.zname=II.Variables(4).Attributes.Value
    else
        INFO.zname=' ';
    end


%     % peruse structure
%     for i=1:nvars
%         natts = numel(II.Variables(i))
%         for j=1:natts
%             II.Variables(i).Attributes(j)
%             %nvals=numel(II.Variables(i).Attributes(j))
%             %for k=1:nvals
%                %nam=II.Variables(i).Attributes(j).Name
%                %val=II.Variables(i).Attributes(j).Value(k)
%             %end
%         end
%     end
    
    fprintf(1,'Attributes.Name\n');
    II.Attributes.Name
    fprintf(1,'Attributes.Value\n');
    II.Attributes(1).Value
    II.Attributes(2).Value

else
    INFO = netcdf.inq(ncid);
    switch nvars
        case 3                       % (v4) GMT netCDF grid file
            x=netcdf.getVar(ncid,0)';
            y=netcdf.getVar(ncid,1)';
            z=netcdf.getVar(ncid,2)';
        case 4                       % v. 5 GMT netCDF grid file
            x=netcdf.getVar(ncid,0)';
            y=netcdf.getVar(ncid,1)';
            z=netcdf.getVar(ncid,2)';
            grid_mapping=netcdf.getVar(ncid,3); %
        case 6                       % old (v3) GMT netCDF grid file
            [dimname, dimlen] = netcdf.inqDim(ncid,1);
            if (dimname=='xysize')            % make sure it really is v3 netCDF
                xrange=netcdf.getVar(ncid,0)';
                yrange=netcdf.getVar(ncid,1)';
                z=netcdf.getVar(ncid,5);
                dim=netcdf.getVar(ncid,4)';
                %             Error using netcdflib
                %The NetCDF library encountered an error during execution of 'inqAtt'
                %function - 'Attribute not found (NC_ENOTATT)'.

                pixel=netcdf.getAtt(ncid,5,'node_offset')
                if pixel                         % pixel node registered
                    dx=diff(xrange)/double(dim(1)); % convert int to double for division
                    dy=diff(yrange)/double(dim(2));
                    x=xrange(1)+dx/2:dx:xrange(2)-dx/2; % convert to gridline registered
                    y=yrange(1)+dy/2:dy:yrange(2)-dy/2;
                else                              % gridline registered
                    dx=diff(xrange)/double(dim(1)-1); % convert int to double for division
                    dy=diff(yrange)/double(dim(2)-1);
                    x=xrange(1):dx:xrange(2);
                    y=yrange(1):dy:yrange(2);
                end
                z=flipud(reshape(z,dim(1),dim(2))');
            else
                error('Apparently not a GMT netCDF grid');
            end
        otherwise
            error('Wrong number of variables %d in netCDF file.\n',nvar);
    end
    netcdf.close(ncid);

    if pixel                         % pixel node registered
        dx=diff(xrange)/double(nx); % convert int to double for division
        dy=diff(yrange)/double(ny);
    else                              % gridline registered
        dx=diff(xrange)/double(nx-1); % convert int to double for division
        dy=diff(yrange)/double(ny-1);
    end

    INFO.title = title;
    INFO.conventions = conv;
    INFO.gmtversion = vers;
    INFO.command = command;
    INFO.description = desc;
    INFO.ispixelreg = pixel;
    INFO.dx = dx;
    INFO.dy = dy;
    INFO.xmin = xrange(1);
    INFO.xmax = xrange(2);
    INFO.ymin = yrange(1);
    INFO.ymax = yrange(2);
    INFO.zmin = zrange(1);
    INFO.zmax = zrange(2);
    INFO.nx   = nx;
    INFO.ny   = ny;
    INFO.xname = xname;
    INFO.yname = yname;
    INFO.zname = zname;

    % disp(['Title: ' title]);
    % disp(['Conventions: ' conv]);
    % disp(['GMT version: ' vers]);
    % disp(['Command: ' command]);
    % disp(['Remark: ' desc]);
    % if pixel,
    %   disp('Pixel node registration used');
    % else
    %   disp('Gridline node registration used');
    % end
    % disp(['x_min: ' num2str(xrange(1)) ' x_max: ' num2str(xrange(2)) ...
    %                    ' name: ' xname])
    % disp(['y_min: ' num2str(yrange(1)) ' y_max: ' num2str(yrange(2)) ...
    %                    ' name: ' yname])
    % disp(['z_min: ' num2str(zrange(1)) ' z_max: ' num2str(zrange(2)) ...
    %                    ' name: ' zname])
    %
    % switch nargout
    %   case 1,
    %    % need to make sure that zrange and pixel are cast into double precision
    %    d=[xrange(1) xrange(2) yrange(1) yrange(2) double(zrange(1)) double(zrange(2)) double(pixel) dx dy];
    % end

end

return

function val = getatt_clean(ncid,attname)
%
% Call netcdf.getAtt when not sure whether global attribute 'attname' exists.
% Trap error and return empty vector in that case.
%
eval(['val = netcdf.getAtt(ncid,netcdf.getConstant(''NC_GLOBAL''),''' attname ''');'],'val=[];');

return


