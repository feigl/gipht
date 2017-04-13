function plot_displacement(DISPLACEMENT,UNITV,BBOX)%% given results of model, plot resutls
% Kurt Feigl 20170116
%initialize
nf = 0;




%% plot the vertical displacement in model coordinates
nf=nf+1;h(nf)=figure;hold on;
imagesc(reshape(DISPLACEMENT.uv,DISPLACEMENT.nrows,DISPLACEMENT.ncols));
colormap(jet);
xlabel('index of Xmodel');
ylabel('index of Ymodel');
axis tight;
title('Vertical displacement [m]');
colorbar;
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%% map of observation points in UTM coordinates
nf=nf+1;h(nf)=figure;hold on;axis equal; axis xy;
plot(DISPLACEMENT.e/1.e3,DISPLACEMENT.n/1.e3,'k.');
xlabel('Easting [km]');
ylabel('Northing [km]');
title('DISPLACEMENT.e and DISPLACEMENT.n');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%% interpolate to regular grid in UTM easting and northing
% griddata Interpolates scattered data - generally to produce gridded data
%     Vq = griddata(X,Y,V,Xq,Yq)
% map coordinates
% pixel dimensions in easting and northing
de = 50; % meters
dn = 50; % meters
%[M.e,M.n] = meshgrid([min(DISPLACEMENT.e):de:max(DISPLACEMENT.e)],[min(DISPLACEMENT.n):dn:max(DISPLACEMENT.n)]);
% [M.e,M.n] = meshgrid([326.6e3:de:329.1e3],[4405.3e3:dn:4408.7e3]);
[M.e,M.n] = meshgrid([BBOX.w:de:BBOX.e],[BBOX.s:dn:BBOX.n]);
M.ue = griddata(DISPLACEMENT.e,DISPLACEMENT.n,DISPLACEMENT.ue,M.e,M.n);
M.un = griddata(DISPLACEMENT.e,DISPLACEMENT.n,DISPLACEMENT.un,M.e,M.n);
M.uv = griddata(DISPLACEMENT.e,DISPLACEMENT.n,DISPLACEMENT.uv,M.e,M.n);

% UNIT VECTOR FROM TARGET TO SATELLITE
% 12 valid records from orbit file /data/bradys/TSX/orbits/headers/6984.orb
% taken from /data/bradys/TSX/raw/20150714_orb.txt
% 0.556398 , -0.101079 , 0.824745
% unitv_east   =  0.556398
% unitv_north  = -0.101079
% unitv_up     =  0.824745

% unitv_east   =  0
% unitv_north  =  0
% unitv_up     =  0.824745

%% range change
%M.ur = -1* (M.ue * unitv_east + M.un * unitv_north + M.uv * unitv_up);
M.ur = -1* (M.ue * UNITV.e + M.un * UNITV.n + M.uv * UNITV.u);

%% clip values
% clip=nan(2,1);
% quantile(M.ur,0.01)
% clip(1) = max(quantile(M.ur,0.01));
% clip(2) = min(quantile(M.ur,0.99));
% ilo = find(M.ur < clip(1));M.ur(ilo) = clip(1);
% ihi = find(M.ur > clip(2));M.ur(ihi) = clip(2);


% ibad = find(abs(M.ur) > 0.02);
% M.ur(ibad) = nan;

%% histogram
nf=nf+1;h(nf)=figure;hold on;
histogram(colvec(M.ur*1.e3),100);
xlabel('range change in mm');
ylabel('number of occurences');
title('histogram of range change');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%% profiles
% find max range change
[iprof,jprof] = find(abs(M.ur - max(max(M.ur))) <= eps);

% profile through max range change on E-W line
nf=nf+1;h(nf)=figure;hold on;
plot(M.e(iprof,:)/1.e3,M.ur(iprof,:)*1.e3,'r+-');
xlabel('Easting [km]');
ylabel('range change in mm');
title(sprintf('Profile on E-W line through Y = %10.3f km',M.n(iprof,jprof)/1.e3));
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

% profile through max range change on N-S line
nf=nf+1;h(nf)=figure;hold on;
plot(M.n(:,jprof)/1.e3,M.ur(:,jprof)*1.e3,'r+-');
xlabel('Northing [km]');
ylabel('range change in mm');
title(sprintf('Profile on N-S line through X = %10.3f km',M.e(iprof,jprof)/1.e3));
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));




%% color
nf=nf+1;h(nf)=figure;hold on;axis equal; axis xy; axis tight;
colormap(jet)
imagesc(linspace(min(min(M.e)),max(max(M.e)),5)/1.e3 ...
    ,linspace(min(min(M.n)),max(max(M.n)),5)/1.e3 ...
    ,M.ur*1.e3);
[C,h] = contour(M.e/1e3,M.n/1e3,M.ur*1e3,'w-');
clabel(C,h,'LabelSpacing',72,'Color','w','FontWeight','bold')
colorbar;
xlabel('Easting [km]');
ylabel('Northing [km]');
title('range change in mm');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%% quiver
nf=nf+1;h(nf)=figure;hold on;axis equal; axis xy; axis tight;
quiver(M.e/1.e3,M.n/1.e3,M.ue*1.e3,M.un*1e3);
xlabel('Easting [km]');
ylabel('Northing [km]');
title('Horizontal displacement vectors in mm');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));


return










