function varargout =grdread4(file)
%function [x,y,z,d]=grdread4(file)
%GRDREAD3  Load a GMT grdfile (netcdf format)
%
% Uses NetCDF libraries to load a GMT grid file.
% Duplicates (some) functionality of the program grdread (which requires
% compilation as a mexfile-based function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% Z=GRDREAD2('filename.grd') will return the data as a matrix in Z
%
% [X,Y,Z]=GRDREAD2('filename.grd') will also return X and Y vectors
% suitable for use in Matlab commands such as IMAGE or CONTOUR.
% e.g., imagesc(X,Y,Z); axis xy
%
% Although both gridline and pixel registered grids can be read,
% pixel registration will be converted to gridline registration
% for the x- and y-vectors.
%
% See also GRDWRITE3, GRDINFO3
%
% CAUTION: This program currently does little error checking and makes
% some assumptions about the content and structure of NetCDF files that
% may not always be valid.  It is tested with COARDS-compliant NetCDF
% grdfiles, the standard format in GMT 4 and later, as well as GMT v3
% NetCDF formats.  It will not work with any binary grid file formats.
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
% Time-stamp: <Wed Jan  6 16:37:45 EST 2010>
%
% Version 1.1.1, 6-Jan-2010
% released with minor changes in documentation along with grdwrite2 and grdinfo2
% Version 1.1, 3-Dec-2009
% support for GMT v3 grids added
% Version 1.0, 29-Oct-2009
% first posted on MATLAB Central
% attempt to adapt to GMT 5 Kurt Feigl 20160813
% 2021/03/10 Kurt Feigl handle error about attributes
% 2022/08/23 Kurt Feigl use netcdf functions 

if nargin < 1
    help(mfilename);
    return,
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
    %% new style
    INFO = ncinfo(file);
    nvars=numel(INFO.Variables);

    if INFO.Variables(1).Attributes(4).Name ==  'node_offset'
       pixel=INFO.Variables(1).Attributes(4).Value;
    else
       fprintf(1,'Assuming grid-node registration.\n');
       pixel=0;
    end
    if nvars == 6
        for i=1:nvars
            name1=INFO.Variables(i).Name
            val=ncread(file,name1);
            switch name1
                case {'z'}
                    z=val;
                case {'x_range'}
                    x_range=val;
                case {'y_range'}
                    y_range=val;
                case {'z_range'}
                    z_range=val;
                case {'spacing'}
                    spacing=val;
                    dx=spacing(1);
                    dy=spacing(2);
                case {'dimension'}
                    dim=val;
                    nx = dim(1);
                    ny = dim(2);
                otherwise
                    warning(sprintf('Unknown variable %s',name1));
            end
        end
    else
        error(sprintf('incorrect dimension %d',6));
    end
   
    %% check sanity
    nz = numel(z);
    if  nz ~= nx * ny
        error(sprintf('nx = %d ny = %d does not equal nz = %d',nx,ny,nz));
    end
    if pixel == 1
        if    abs(x_range(1)+(nx-1)*dx - x_range(2)) < dx/2. ...
           && abs(y_range(1)+(ny-1)*dy - y_range(2)) < dy/2.     
            fprintf(1,'spacing for pixel-centered registration OK\n')            
        else
            error(sprintf('spacing issue for pixel-centered registration\n'));
        end
    else
        if    abs(x_range(1)+ nx*dx - x_range(2)) < dx/2. ...
           && abs(y_range(1)+ ny*dy - y_range(2)) < dy/2.
            fprintf(1,'spacing for grid-node registration OK\n')
        else
            error(sprintf('spacing issue grid-node registration\n'));
        end
    end
    if pixel                         % pixel node registered
        ddx=diff(x_range)/double(dim(1)); % convert int to double for division
        ddy=diff(y_range)/double(dim(2));
        x=x_range(1)+ddx/2:ddx:x_range(2)-ddx/2; % convert to gridline registered
        y=y_range(1)+ddy/2:ddy:y_range(2)-ddy/2;
    else                              % gridline registered
        ddx=diff(x_range)/double(dim(1)-1); % convert int to double for division
        ddy=diff(y_range)/double(dim(2)-1);
        x=x_range(1):ddx:x_range(2);
        y=y_range(1):ddy:y_range(2);
    end
    z=flipud(reshape(z,dim(1),dim(2))');
    whos
else
    %% old-style reading for backward compatibility
    ncid = netcdf.open(file, 'NC_NOWRITE');
    if isempty(ncid)
        return;
    end

    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    nvars
    ngatts
    ndims
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
                x_range=netcdf.getVar(ncid,0)';
                y_range=netcdf.getVar(ncid,1)';
                z=netcdf.getVar(ncid,5);
                dim=netcdf.getVar(ncid,4)';
                %             Error using netcdflib
                %The NetCDF library encountered an error during execution of 'inqAtt'
                %function - 'Attribute not found (NC_ENOTATT)'.

                pixel=netcdf.getAtt(ncid,5,'node_offset')
                if pixel                         % pixel node registered
                    dx=diff(x_range)/double(dim(1)); % convert int to double for division
                    dy=diff(y_range)/double(dim(2));
                    x=x_range(1)+dx/2:dx:x_range(2)-dx/2; % convert to gridline registered
                    y=y_range(1)+dy/2:dy:y_range(2)-dy/2;
                else                              % gridline registered
                    dx=diff(x_range)/double(dim(1)-1); % convert int to double for division
                    dy=diff(y_range)/double(dim(2)-1);
                    x=x_range(1):dx:x_range(2);
                    y=y_range(1):dy:y_range(2);
                end
                z=flipud(reshape(z,dim(1),dim(2))');
            else
                error('Apparently not a GMT netCDF grid');
            end
        otherwise
            error('Wrong number of variables %d in netCDF file.\n',nvar);
    end
    netcdf.close(ncid);
end

%% return variables
switch nargout
    case 1
        varargout{1}=z;
    case {3,4}
        varargout{1}=x;
        varargout{2}=y;
        varargout{3}=z;
        if nargout == 4 && nvars == 4
            varargout{4}=grid_mapping;
        end
    otherwise
        error('grdread4: Incorrect number of output arguments!');
end
return

