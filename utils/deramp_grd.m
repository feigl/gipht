function [fgridout,mest,msig,mse] = deramp_grd(xvec,yvec,fgridin)
% given an rectangular array, estimate 3 parameters and remove a ramp
% 20190220 Kurt Feigl

[ny,nx] = size(fgridin);

if ny == numel(yvec) || nx == numel(xvec)


    
% find coordinates of all points with respect to center
[xgrid,ygrid] = meshgrid(xvec,yvec);
xgrid = colvec(xgrid - nanmean(xvec));
ygrid = colvec(ygrid - nanmean(yvec));

fdata = colvec(fgridin);
iok = find(isfinite(fdata)==1);

ndat = numel(iok);



% design matrix
GG = zeros(ndat,3);
for i=1:ndat
    GG(i,1) = 1.0; % additive constant
    GG(i,2) = xgrid(iok(i)); % X component of gradient
    GG(i,3) = ygrid(iok(i)); % Z component of gradient
end

%% data weights
% W = ones(ndat,1)/std(fdata(iok));
%% perform least squares estimate 
[ mest, msig, mse, S] = lscov(GG,fdata(iok));

fprintf(1,'sqrt(Mean Squared Error) = %12.4g\n', sqrt(mse));
fprintf(1,'Parameter         Estimate    Uncertainty\n');
fprintf(1,'Additive_Constant %10.4g  +/- %10.4g\n',mest(1),msig(1)); 
fprintf(1,'X-gradient        %10.4g  +/- %10.4g\n',mest(2),msig(2)); 
fprintf(1,'Y-gradient        %10.4g  +/- %10.4g\n',mest(3),msig(3)); 

fmodel = mest(1) + mest(2)*xgrid + mest(3)*ygrid;

fresid = fdata - fmodel;
figure
histogram(fresid);
close

fgridout = reshape(fmodel,ny,nx);

else
    nx
    ny
    error('miscount');
end

return
end

