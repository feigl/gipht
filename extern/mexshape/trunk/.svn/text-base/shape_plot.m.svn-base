function h = shape_plot ( shape_struct )
% SHAPE_PLOT:  quick and dirty plotting routine for data read by shape_in
%
% USAGE:  h = shape_plot ( shapestruct, shape_type );
%
% PARAMETERS:
% Input:
%     shape_struct:
%         structure of shapes read by shape_in
%     shape_type:
%         string determining the type of shapefile, also returned by shape_in
% Output:
%     h:
%         vector of handles, one for each element of shape_struct
%


%
% If empty, there is nothing to do
if length(shape_struct) < 1
	return
end

oldNextPlot = get ( gca, 'NextPlot' );

n = length(shape_struct.Shape);
for j = 1:n
	hold on
	switch ( shape_struct.Type )
		case {'Arc', 'MultiPoint', 'Point'}
			h(j,1) = plot ( shape_struct.Shape(j).x, shape_struct.Shape(j).y );
		case 'Polygon'
			ind = find(isnan(shape_struct.Shape(j).x));
			part_start = [1; ind(1:end-1)+1];
			part_stop = ind-1;
			for k = 1:length(ind)
				index_range = part_start(k):part_stop(k);
				h(j,k) = patch ( shape_struct.Shape(j).x(index_range), shape_struct.Shape(j).y(index_range), 'r' );
			end
	end

end

set ( gca, 'NextPlot', oldNextPlot );
