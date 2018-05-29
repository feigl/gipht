function [s, t] = mex_shape ( shapefile )
% MEX_READSHAPEFILE:  reads in an ESRI shapefile and attributes
%
% USAGE:  [s, t] = mex_shape ( shapefile );
%
% PARAMETERS:
% Input:
%    shapefile:
%        path of shapefile.  This can end in the .shp suffix or possibly not.
% Output:
%    s:
%        structure array of shapes.  Each structure has as fields each of the
%        attributes stored in the corresponding record in the DBF file, as well
%        as two fields "mx_data" and "my_data", which contain the xy coordinates.
%        If the DBF file has fields called "mx_data" or "my_data", we're screwed.
%    t:  
%        type of shapes read in.  This can be 'Polygon', 'Point', or 'Arc'
%
% This mex routine should not be called directly.  Use the convenience function 
% "shape_in" instead.
%
% AUTHOR:
%    johnevans@acm.org

