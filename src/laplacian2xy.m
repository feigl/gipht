function LL = laplacian2xy(x,y)
%% Form 2-D finite difference for two-dimensional Laplacian operator
% input: nx, ny number of points in x and y
% output: LL sparse, square matrix with dimensions (nx*ny) by (nx*ny)
% 20170515 Kurt Feigl

nx = numel(x);
ny = numel(y);
if nx == ny   
    ncells = nx
else
    error
end
% number of nonzero entries in sparse matrix
nxy = nx*ny
% number of nonzero entries in sparse matrix
nzmax = 5*(nx-2)*(ny-2)


%% make convex hull in 2 dimensions
DelTriang = delaunayTriangulation(x,y);
khull = convexHull(DelTriang);
hullx = x(khull);
hully = y(khull);

%% find the points on the edge
[iin,ion] = inpolygon(x,y,hullx,hully)

non = numel(find(ion==1))

% number of nonzero rows
nzrows = ncells - non
% note order!
LL=zeros(nzrows,nxy);
size(LL)


figure;hold on;
plot(x,y,'ro');
plot(hullx,hully,'b-');
plot(x(ion),y(ion),'k*');

kount = 0;
k=0;
% for irow = 1:nxy
%     [j,i] = ind2sub([ncells,ncells],irow);
for ii=1:ncells
    for jj=1:ncells
        k = k+1;
        [j,i] = ind2sub([ncells,ncells],k);
        if i>1 && i<ncells && j>1 && j<ncells
            if ion(j-1)==0 && ion(j)==0 && ion(j+1) && ion(i-1)==0 && ion(i)==0 && ion(i+1)==0
                kount = kount+1;
                %         fprintf(1,'%3dth nonzero row at interior (j,i) %3d %3d\n',kount,j,i);
                icol = sub2ind([ncells,ncells],j  , i  );
                LL(kount,icol) = -1;   % central voxel
                icol = sub2ind([ncells,ncells],j  , i+1); LL(kount,icol) =  1/4; % right neighbor in X
                icol = sub2ind([ncells,ncells],j  , i-1); LL(kount,icol) =  1/4; % left  neighbor in X
                icol = sub2ind([ncells,ncells],j+1, i  ); LL(kount,icol) =  1/4; % right neighbor in Y
                icol = sub2ind([ncells,ncells],j-1, i  ); LL(kount,icol) =  1/4; % left  neighbor in Y
            end
        end
    end
end
figure;
spy(LL)

if kount ~= nzrows
    kount
    nzrows
    nzmax
    error('miscount');
end

return

end


