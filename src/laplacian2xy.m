function LL = laplacian2xy(x,y)
%% Form 2-D finite difference for two-dimensional Laplacian operator
% input: nx, ny number of points in x and y
% output: LL sparse, square matrix with dimensions (nx*ny) by (nx*ny)
% 20170515 Kurt Feigl

nx = numel(x);
ny = numel(y);
if nx == ny   
    ncells = nx;
else
    error
end
% number of nonzero entries in sparse matrix
LL=zeros(ncells,ncells);

%% make convex hull in 2 dimensions
DelTriang = delaunayTriangulation(x,y);
khull = convexHull(DelTriang);
hullx = x(khull);
hully = y(khull);

%% find the points on the edge
[iin,ion] = find(inpolygon(x,y,hullx,hully));

kount = 0;
for irow = 1:ncells
    [j,i] = ind2sub([ny,nx],irow);
    if ismember(irow,ion) == 0
        kount = kount+1;
        fprintf(1,'%3dth nonzero row at interior (j,i) %3d %3d\n',kount,j,i);
        
        icol = sub2ind([ny,nx],j  , i  ); LL(irow,icol) = -1;   % central voxel
        icol = sub2ind([ny,nx],j  , i+1); LL(irow,icol) =  1/4; % right neighbor in X
        icol = sub2ind([ny,nx],j  , i-1); LL(irow,icol) =  1/4; % left  neighbor in X
        icol = sub2ind([ny,nx],j+1, i  ); LL(irow,icol) =  1/4; % right neighbor in Y
        icol = sub2ind([ny,nx],j-1, i  ); LL(irow,icol) =  1/4; % left  neighbor in Y
    elseif i == 1 || i == nx || j == 1 || j == ny
        kount = kount+1;
        fprintf(1,'%3dth nonzero row at edge     (j,i) %3d %3d\n',kount,j,i);        
        icol = sub2ind([ny,nx],j  , i  ); LL(irow,icol) = 1;   % central voxel
    end
end

return

end


