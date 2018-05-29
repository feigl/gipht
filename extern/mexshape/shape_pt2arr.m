function [xdata, ydata] = shape_pt2arr ( shape_struct )
% SHAPE_PT2ARR:  extracts points from a shapefile into arrays
% 
% USAGE:  [x, y] = shape_pt2arr ( shape_struct );
%
% PARAMETERS:
% Input:
%    shape_struct:
%        structure returned by mex_shape
% Output:
%    x, y:  double precision arrays consisting of the point values 
%
% AUTHOR:  
% John Evans, johnevans@acm.org


n = length(shape_struct.Shape);
xdata = NaN * ones(n,1);
ydata = NaN * ones(n,1);
for j = 1:n
	xdata(j,1) = shape_struct.Shape(j).x;
	ydata(j,1) = shape_struct.Shape(j).y;
end


