function Qiev = edges_to_incidence(edges)
% function [Qiev, trees] = edges_to_incidence(edges)
% inputs:
%    edges == list of pairs of indices to vertices, one row per pair
% outputs: 
%     Qiev  == edge-vertex incidence matrix
%
% simplified from find_trees by Kurt L. Feigl, Elena C. Reinisch, UW-Madison
% 
% Reference:
% Reinisch, E.C., Cardiff, M. & Feigl, K.L. Graph theory for analyzing
% pair-wise data: application to geophysical model parameters estimated
% from interferometric synthetic aperture radar data at Okmok volcano,
% Alaska. J Geod 91, 9?24 (2017).
% https://doi.org/10.1007/s00190-016-0934-5
%
% 201200326


narginchk(1, 1);
nargoutchk(1, 2);

% verbose ?
dispflag = 0;


if dispflag == 1
    fprintf(1,'%s begins ...\n',mfilename);
end

% initialize values to return
Qiev = NaN;
ivertices = NaN;
tepochs = NaN;

[nedges, ncols] = size(edges);
if ncols ~= 2
    ncols
    error(sprintf('input matrix edges must have only 2 columns. It has %d\n',ncols));
end

if dispflag == 1
    disp 'number of pairs'
    nedges
end

% 
%% integer indices to vertices
i1 = edges(:,1);
i2 = edges(:,2);


[ivertices,iuniq,juniq] = unique([i1 i2]);
 
% find indices
ibegins = zeros(nedges,1);
ifinish = zeros(nedges,1);
for i=1:nedges
    ibegins(i) = find(ivertices == i1(i)); % index of starting vertex
    ifinish(i) = find(ivertices == i2(i)); % index of finishing vertex
end
    
% Set number of vertices
nvertices = length(ivertices);

% Display values if desired
if dispflag == 1
    disp 'number of distinct epochs', nvertices
    disp 'find differencing operator Q such that Q*epochs = pairs ...'
end

% Initialize edge-vertex incidence matrix Qiev
Qiev = zeros(nedges,nvertices);
for p = 1:nedges
   for j = 1:nvertices
      if i1(p) == ivertices(j)
         Qiev(p,j) = -1; % assign -1 if j_th epoch is the master epoch of pair p
      elseif i2(p) == ivertices(j) 
         Qiev(p,j) = 1; % assign +1 if j_th epoch is the slave epoch of pair p
      end
   end
end

if dispflag == 1
   Qiev
end

return;

end



