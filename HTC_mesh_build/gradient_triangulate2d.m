function DEL = gradient_triangulate2d(x,y)
%% Form 2-D finite difference for gradient operator in 2 dimensions by triangulation
% input: x, y - column vectors of coordinates
% output: DEL square matrix with dimensions (nx*ny) by (nx*ny)
% 20170515 Kurt Feigl

nx = numel(x);
ny = numel(y);
if nx == ny   
    ncells = nx
else
    error
end
% number of nonzero entries in sparse matrix
% nxy = nx*ny
% % number of nonzero entries in sparse matrix
% nzmax = 5*(nx-2)*(ny-2)

%% compute Delaunay Triangulation
DT = delaunayTriangulation(x,y)
whos DT
TR = triangulation(DT.ConnectivityList,x,y)
whos TR
[ntriangles,ncols] = size(TR)

%% find list of edges
kedges = edges(TR);
[nedges,ncols] = size(kedges);

%% make convex hull in 2 dimensions
khull = convexHull(DT);
hullx = x(khull);
hully = y(khull);

%% find the points on the edge
[iin,ion] = inpolygon(x,y,hullx,hully);
non = numel(find(ion==1))



figure;hold on;
plot(x,y,'ro');
plot(hullx,hully,'b-');
plot(x(ion),y(ion),'k*');
for i=1:nedges
   plot([x(kedges(i,1));x(kedges(i,2))],[y(kedges(i,1));y(kedges(i,2))],'g--');
end
% for i=1:ntriangles
%     plot(x(DT(i,:)),y(DT(i,:)),'r:');
% end
axis xy
axis equal


% number of nonzero rows
nzrows = nedges;
% note order!
DEL=zeros(nzrows,ncells);
size(DEL)

kount = 0;
k=0;
% for irow = 1:nxy
%     [j,i] = ind2sub([ncells,ncells],irow);
for jj=1:nedges
    i1=kedges(jj,1);
    i2=kedges(jj,2);
    if ion(i1) == 0 && ion(i2) == 0
        x1=x(i1);
        x2=x(i2);
        y1=y(i1);
        y2=y(i2);
        dist=hypot(x2-x1,y2-y1);       
        DEL(jj,i1) = -1./dist;
        DEL(jj,i2) = +1./dist;
    elseif ion(i1) == 1
        DEL(jj,i1) = 1;
    elseif ion(i2) == 1
        DEL(jj,i1) = 1;
    end
end
DEL=sparse(DEL);
figure;
spy(DEL)

% if kount ~= nzrows
%     kount
%     nzrows
%     nzmax
%     error('miscount');
% end

return

end



