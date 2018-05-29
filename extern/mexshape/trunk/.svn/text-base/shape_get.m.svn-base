function [shape_data, xpoint, ypoint] = shape_get ( shapefile )
% SHAPE_IN:  reads an ESRI shapefile and associated attributes.
%
% So far, this routine is limited to use with point, polyline, and
% polygon shapefiles.  Other types (and there are lots of others)
% are not handled and will cause an exception to be thrown.
%
% USAGE:
%     [shape_data] = shape_get ( shapefile );
%
%     If the shape type is "Point", then the points can optionally
%     be put into matlab arrays for easier access with a call like
%
%     [shape_data, x, y] = shape_get ( shapefile );
%
% PARAMETERS:
% Input:
%    shapefile:
%        A shapefile usually consists of three files.  As stated by Frank Warmerdam,
%        the author of shapelib, the purpose of each of the three files is
%
%           XXX.shp - holds the actual data vertices
%           XXX.shx - holds index data pointing to the structures in the .shp file
%           XXX.dbf - holds the attributes in xBase (dBase) format.
%
%        The argument given for "shapefile" here can be either the full path to the
%        XXX.shp file, or one can do away with the ".shp" part altogether.
%
% Output:
%    shape_data:
%        structure of shape data.  Each structure has three fields,
%            Shape:
%                Array of structures with "x" and "y" fields that
%                correspond to the xy coordinates in the shapefile.
%            Attribute:
%                Array of structure of attributes.  Each element
%                corresponds to one with the same index in the 
%                "Shape" field, and each such element has as its
%                own fields the attributes which are stored in the
%                DBF file.
%            Type:
%                String with four possible values, "Arc", "MultiPoint", 
%                "Point", or "Polygon".
%    x, y:
%        Optional.  If given and if the data is 'Point', then all 
%        of the point values are extracted to single arrays. 
%
% AUTHOR: 
% John Evans, john.g.evans.ne@gmail.com


xpoint = [];
ypoint = [];

if nargout < 1
	fprintf ( 2, '%s:  at least one output argument is required.\n', mfilename );
	help shape_get;
	return
end

if ~exist ( shapefile, 'file' )
	fprintf ( 2, '%s:  The file ''%s'' does not exist.\n', mfilename, shapefile );
	return
end

[shape_data] = mexshape ( shapefile );
switch shape_data.Type

	case 'Point'

		%
		% If the calling sequence was [x, y, .. ] = shape_get ( .. );
		% then extract the point data into single arrays.
		if ( nargout >= 1 )
			[xpoint, ypoint] = shape_pt2arr ( shape_data );
		end


	case { 'Arc', 'MultiPoint', 'Polygon' }


	otherwise
		msg = sprintf ( '%s:  unhandled shape type %s.\n', mfilename, shape_type );
		error ( msg );
end

