function LL = laplacian2(ny,nx)
%% Form 2-D finite difference for two-dimensional Laplacian operator
% input: nx, ny number of points in x and y
% output: LL sparse, square matrix with dimensions (nx*ny) by (nx*ny)
% 20170515 Kurt Feigl

nxy = nx*ny;
% number of nonzero entries in sparse matrix
nzmax = 5*((nx-2)+(ny-2))
%LL = spalloc(nxyz,nxyz,nzmax);
% note order!
LL=zeros(nxy,nxy);
kount = 0;
for irow = 1:nxy
    [j,i] = ind2sub([ny,nx],irow);
    if i>1 && i<nx && j>1 && j<ny 
        kount = kount+1;
        fprintf(1,'%3dth nonzero row at (j,i) %3d %3d\n',kount,j,i);
        icol = sub2ind([ny,nx],j  , i  ); LL(irow,icol) = -1;   % central voxel
        icol = sub2ind([ny,nx],j  , i+1); LL(irow,icol) =  1/4; % right neighbor in X
        icol = sub2ind([ny,nx],j  , i-1); LL(irow,icol) =  1/4; % left  neighbor in X
        icol = sub2ind([ny,nx],j+1, i  ); LL(irow,icol) =  1/4; % right neighbor in Y
        icol = sub2ind([ny,nx],j-1, i  ); LL(irow,icol) =  1/4; % left  neighbor in Y
    end
end

return

end


