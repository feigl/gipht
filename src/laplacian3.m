function LL = laplacian3(nx,ny,nz)
%% Form 3-D finite difference for three-dimensional Laplacian operator
% input: nx, ny, nz, number of points in x, y, z
% output: LL sparse, square matrix with dimensions (nx*ny*nz) by (nx*ny*nz)
% 20170512 Kurt Feigl

nxyz = nx*ny*nz;
% number of nonzero entries in sparse matrix
nzmax = 7*((nx-2)+(ny-2)+(nz-2))
%LL = spalloc(nxyz,nxyz,nzmax);
% note order!
LL=zeros(nxyz,nxyz);
kount = 0;
for irow = 1:nxyz
    [k,j,i] = ind2sub([nz,ny,nx],irow);
    if i>1 && i<nx && j>1 && j<ny && k>1 && k<nz
        kount = kount+1;
        fprintf(1,'%3dth nonzero row at (k,j,i) %3d %3d %3d\n',kount,k,j,i);
        icol = sub2ind([nz,ny,nx],k , j  , i  ); LL(irow,icol) = -1;   % central voxel
        icol = sub2ind([nz,ny,nx],k , j  , i+1); LL(irow,icol) =  1/6; % right neighbor in X
        icol = sub2ind([nz,ny,nx],k , j  , i-1); LL(irow,icol) =  1/6; % left  neighbor in X
        icol = sub2ind([nz,ny,nx],k , j+1, i  ); LL(irow,icol) =  1/6; % right neighbor in Y
        icol = sub2ind([nz,ny,nx],k , j-1, i  ); LL(irow,icol) =  1/6; % left  neighbor in Y
        icol = sub2ind([nz,ny,nx],k+1,j  , i  ); LL(irow,icol) =  1/6; % right neighbor in Z
        icol = sub2ind([nz,ny,nx],k-1,j  , i  ); LL(irow,icol) =  1/6; % left  neighbor in Z
    end
end

return

end


