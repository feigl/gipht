function modeledDisplacements = fun_fitmogi3(PST,DST,TST)
% PST parameter structure
% DST data structure
% fitting function
%persistent kount

narginchk(3,3);

% %fprintf(1,'Entering %s\n',mfilename);
% mparams = numel(PST.p1);
% for i=1:mparams
% fprintf(1,'%10.4g ',PST.p1(i));
% end
% fprintf(1,'\n');


%% calculate modeled values at locations of data points
nsites = numel(DST.sites);
volgeom(1) = PST.p1(1);          % E coordinate of source in m
volgeom(2) = PST.p1(2);          % N coordinate of source in m
volgeom(3) = PST.p1(3);          % Depth of source in m
volgeom(4) = PST.p1(4);          % Volume change of source in m^3
nu         = PST.p1(5);          % Poisson ratio, dimensionless 
umogi = mogi(volgeom, [DST.x';DST.y'], nu);

ndata = 3*nsites;

modeledDisplacements = nan(ndata,1);
idata=0;
for isite = 1:nsites
    idata=idata+1;modeledDisplacements(idata) = umogi(1,isite);
    idata=idata+1;modeledDisplacements(idata) = umogi(2,isite);
    idata=idata+1;modeledDisplacements(idata) = umogi(3,isite);
end


if idata ~= ndata
    error('miscount')
end

return
end

