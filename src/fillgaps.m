function [image_out] = fillgaps(image_in)
% replace NaNs with bilinear interpolation

[nrows,ncols] = size(image_in)

figure
imagesc(1:ncols,1:nrows,image_in);
colorbar


% make mesh
[X,Y] = meshgrid(1:ncols,1:rows);
size(X)
size(Y)

% create interpolating function
%Finterp = scatteredInterpolant(colvec(X),colvec(Y),colvec(image_in));
%Finterp = scatteredInterpolant(X,Y,image_in,'nearest','none');

% find gaps
[igaps,jgaps] = find(isnan(image_in)==1);

image_out = image_in;
%image_out(ind2sub([nrows,ncols],igaps)) = Finterp(X,Y);
%image_out(igaps) = Finterp(X(igaps),Y(igaps));
for k=1:numel(igaps)
    if igaps(k) > 1 && igaps(k) < nrows && jgaps(k) > 1 && jgaps(k) < ncols
        [ii,jj] = meshgrid(igaps(k)-1:igaps(k)+1, jgaps(k)-1:jgaps(k)+1);
        image_out(igaps(k),jgaps(k)) = nanmean(image_in(ii,jj));
    end
end

figure
imagesc(1:ncols,1:nrows,image_in);
colorbar

return

end

