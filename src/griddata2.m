function zi = griddata2(x,y,z,xi,yi,interp_method,extrap_method)
% replace old GRIDDATA with TriScatteredInterp
%  GRIDDATA Data gridding and surface fitting.
%
%     GRIDDATA is not recommended. Use TriScatteredInterp instead.
%
%     ZI = GRIDDATA(X,Y,Z,XI,YI) fits a surface of the form Z = F(X,Y) to the
%     data in the (usually) nonuniformly-spaced vectors (X,Y,Z). GRIDDATA
%     interpolates this surface at the points specified by (XI,YI) to produce
%     ZI.  The surface always goes through the data points. XI and YI are
%     usually a uniform grid (as produced by MESHGRID) and is where GRIDDATA
%     gets its name.
%
%     XI can be a row vector, in which case it specifies a matrix with
%     constant columns. Similarly, YI can be a column vector and it specifies
%     a matrix with constant rows.
%
%     [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
%     this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%     [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD) where METHOD is one of
%         'linear'    - Triangle-based linear interpolation (default)
%         'cubic'     - Triangle-based cubic interpolation
%         'nearest'   - Nearest neighbor interpolation
%         'v4'        - MATLAB 4 griddata method
%     defines the type of surface fit to the data. The 'cubic' and 'v4'
%     methods produce smooth surfaces while 'linear' and 'nearest' have
%     discontinuities in the first and zero-th derivative respectively.  All
%     the methods except 'v4' are based on a Delaunay triangulation of the
%     data.
%     If METHOD is [], then the default 'linear' method will be used.
%
% Options for interp_method
%             'natural'   Natural neighbor interpolation
%             'linear'    Linear interpolation (default)
%             'nearest'   Nearest neighbor interpolation
%             'nearnat'   Nearest neighbor at edges, natural elsewhere
% Options for extrap_method
%             'natural'   Natural neighbor extrapolation
%             'linear'    Linear interpolation 
%             'none'      Do not extrapolate (default)

% 2022/09/26 Kurt Feigl


if exist('interp_method','var') == 1
    if strcmpi(method,'v4') == 1
        interp_method = 'natural'; % Natural neighbor interpolation
    end
else
    interp_method = 'linear'; % Natural neighbor interpolation
end

if exist('extrap_method','var') == 1
    if strcmpi(method,'v4') == 1
        extrap_method = 'natural'; % Natural neighbor interpolation
    end
else
    extrap_method = 'none'; % default is not to extrapolate
end


    

fprintf(1,'Starting GRIDDATA2 with interpolation method = %s and extrpolation method =%s\n',interp_method,extrap_method);
t0 = tic;

% make a grid
%if nargout > 1
%    [xi, yi] = meshgrid(xi,yi);
%end

[nrowsx, ncolsx] = size(xi);
[nrowsy, ncolsy] = size(yi);

if nrowsx == nrowsy && ncolsx == ncolsy
    nrows = nrowsx;
    ncols = ncolsx;
else
    error(sprintf('ERROR: mismatched sizes %d %d %d %d\n',nrowsx,nrowsy,ncolsx,ncolsy));
    return
end

if strcmpi(interp_method,'nearnat') == 1
    % construct interpolating function
    F1 = TriScatteredInterp(colvec(x),colvec(y),colvec(z),'natural');    
    F2 = TriScatteredInterp(colvec(x),colvec(y),colvec(z),'nearest'); 
    
    % evaluate interpolating functions at all locations
    zi =F1(xi,yi);
    % evaluate interpolating functions at edges
    %nedge = 16;
    %nedge = ceil(min([nrows,ncols])/100.);
    nedge=2;
    i = rowvec((1:nedge)'*ones(1,ncols));             j = rowvec(meshgrid(1:ncols,1:nedge)); k = sub2ind([nrows,ncols],i,j); zi(k) = F2(xi(k),yi(k)); % top edge
    i = rowvec((nrows-nedge+1:nrows)'*ones(1,ncols)); j = rowvec(meshgrid(1:ncols,1:nedge)); k = sub2ind([nrows,ncols],i,j); zi(k) = F2(xi(k),yi(k)); % bottom edge
    j = rowvec((1:nedge)'*ones(1,nrows));             i = rowvec(meshgrid(1:nrows,1:nedge)); k = sub2ind([nrows,ncols],i,j); zi(k) = F2(xi(k),yi(k)); % top edge
    j = rowvec((ncols-nedge+1:ncols)'*ones(1,nrows)); i = rowvec(meshgrid(1:nrows,1:nedge)); k = sub2ind([nrows,ncols],i,j); zi(k) = F2(xi(k),yi(k)); % bottom edge
else 
    % prune 
    iok=find(~ismissing(x));
    iok=intersect(iok,find(~ismissing(y)));
    iok=intersect(iok,find(~ismissing(z)));
    % construct interpolating function
    %F = TriScatteredInterp(colvec(x(iok)),colvec(y(iok)),colvec(z(iok)),interp_method,extrap_method);
    F = scatteredInterpolant(colvec(x(iok)),colvec(y(iok)),colvec(z(iok)),interp_method,extrap_method);
    
    % evaluate interpolating functions at locations
    zi =F(xi,yi);    
end

if nargout == 1
    zi = reshape(zi,size(xi));
end

fprintf(1,'%s finished interpolating after %#10.4f seconds\n',mfilename,toc(t0));
return

