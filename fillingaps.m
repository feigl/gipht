function [image_out,ngood,ngaps] = fillingaps(image_in,interpolation_method,extrapolation_method)
% given an image, replace NaNs with interpolated values
%function [image_out,ngood,ngaps] = fillingaps(image_in,interpolation_method,extrapolation_method)
% examples:
% [imAfilled,ngood,ngaps] = fillingaps(imA,'nearest','none');
% [imAfilled,ngood,ngaps] = fillingaps(imA,'linear','linear');
% [imAfilled,ngood,ngaps] = fillingaps(imA,'natural','linear');
% 20171213 Kurt Feigl

[nrows,ncols] = size(image_in);

% for debugging
% figure
% imagesc(1:ncols,1:nrows,image_in);
% title('image_in');
% xlabel('column index');
% ylabel('row index');
% colorbar


% make mesh
[X,Y] = meshgrid(1:ncols,1:nrows);
% size(X)
% size(Y)
x = colvec(X);
y = colvec(Y);

% find gaps
igaps = find(isnan(colvec(image_in))==1);
ngaps = numel(igaps);
% find good pixels
igood = find(isnan(colvec(image_in))==0);
ngood = numel(igood);


% create interpolating function
Finterp = scatteredInterpolant(x(igood),y(igood),colvec(image_in(igood)),interpolation_method,extrapolation_method);

% copy all the pixel values
image_out = image_in;

% overwrite the gaps
image_out(igaps) = Finterp(x(igaps),y(igaps));

% for debugging
% figure
% imagesc(1:ncols,1:nrows,image_out);
% title('image_out');
% xlabel('column index');
% ylabel('row index');
% colorbar

return

end

