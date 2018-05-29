function [shape_data, shape_type xpoint, ypoint] = shape_in ( shapefile )
% SHAPE_IN:  reads an ESRI shapefile and associated attributes.
%
% So far, this routine is limited to use with point, polyline, and
% polygon shapefiles.  Other types (and there are lots of others)
% are not handled and will cause an exception to be thrown.
%
% This routine is for backwards compatibility purposes.  
%
% USAGE:
%     [shape_data, shapetype] = shape_in ( shapefile );
%
%     If the shape type is "Point", then the points can optionally
%     be put into matlab arrays for easier access with a call like
%
%     [shape_data, shapetype, x, y] = shape_in ( shapefile );
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
%        structure array of shapes.  Each structure has as fields each of the
%        attributes stored in the corresponding record in the DBF file, as well
%        as two fields "mx_data" and "my_data", which contain the xy coordinates.
%        If the DBF file has fields called "mx_data" or "my_data", we're in trouble.
%    shape_type:
%        String with three possible values, "Arc", "Point", or "Polygon".
%
% AUTHOR: 
% John Evans, johnevans@acm.org


xpoint = [];
ypoint = [];

if nargout < 1
	fprintf ( 2, '%s:  at least one output argument is required.\n', mfilename );
	help shape_in;
	return
end

if ~exist ( shapefile, 'file' )
	fprintf ( 2, '%s:  The file ''%s'' does not exist.\n', mfilename, shapefile );
	return
end

[shape_data, shape_att, shape_type, xpoint, ypoint] = shape_get ( shapefile );

%
% Fold the "x" and "y" fields of shape_data info shap_att
for k = 1:length(shape_att)
	shape_att(k).mx_data = shape_data.x;
	shape_att(k).my_data = shape_data.y;
end

shape_data = shape_att;

return
