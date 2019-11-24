function ic = scale_color(ctab,zval,zmin,zmax)
%function ic = scale_color(ctab,zval,zmin,zmax)
% given a color table, and values, return the corresponding index into the
% color table
% 20140610 Kurt Feigl

[ncolors,ndum] = size(ctab);

cuts = linspace(zmin,zmax,ncolors+1);

zval = colvec(zval);

% by default, set index to lowest value in color table
ic = ones(size(zval));


% choose the color
igood = find(isfinite(zval)==1);

for i = 1:numel(igood)
    ztemp = zval(igood(i))
    for j=1:ncolors
        if ztemp >= cuts(j) && ztemp < cuts(j+1)
            ic(i) = j;
        end
    end
end

return



