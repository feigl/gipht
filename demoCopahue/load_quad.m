%% load quad-tree LOS rates 
% De: "Paul R Lundgren (329A)" <Paul.R.Lundgren@jpl.nasa.gov> Para: "MARIA
% CORDOVA" <maria.cordova@sernageomin.cl> Enviados: Martes, 26 de Mayo 2015
% 14:16:25 Asunto: Re: Copahue interferogram
% 
% Hi Loreto,
% 
% Here are three down-sampled time series results (two from COSMO-SkyMed ?
% CSK; and one from RADARSAT-2 ? RSAT2). Column values are: vel_los (m/yr),
% rms_vel, center_lon, center_lat, 4 quad corners lon, 4 quad corners lat,
% ground incidence angle, satellite heading (radar is right looking)
% 
% Let me know if you have any questions.
% 
% Regarding GPS processing, just send an email to her with the RINEX files
% (or an ftp site from which to download the data). If I knew you were
% going to be in the US we could have tried to host you at JPL.
% 
% Paul ____________________________ Paul Lundgren, Ph.D. Jet Propulsion
% Laboratory California Institute of Technology 4800 Oak Grove Drive
% Pasadena, CA 91109, USA tel:  +01 818 354-1795 cell: +01 626 375-3638
% fax: +01 818 354-9476 paul.r.lundgren@jpl.nasa.gov

close all;
clear all;
format long; 
format compact;

qfname = 'Quad_ASC_CSK_2013.75_2014_04_09.out';
Q = load(qfname);
vlos       = Q(:,1);
slos       = Q(:,2);
CenterLon  = Q(:,3);
CenterLat  = Q(:,4);
CornersLon = Q(:,5:8);
CornersLat = Q(:,9:12);
incidence  = Q(:,13);
heading    = Q(:,14);

nquads = numel(vlos)

%% find spacing between quads and build grid in geographic coordinates
dUniqLon = unique(diff(unique(colvec(CornersLon))))
dLon = min(dUniqLon)
dUniqLon / dLon
dUniqLat = unique(diff(unique(colvec(CornersLat))))
dLat = min(dUniqLat)
dUniqLat / dLat

Lonmin = min(colvec(CornersLon));
Lonmax = max(colvec(CornersLon));
Latmin = min(colvec(CornersLat));
Latmax = max(colvec(CornersLat));

[LonI,LatI] = meshgrid([Lonmin:dLon:Lonmax],[Latmin:dLat:Latmax]);

%% grid using nearest neighbor
%vlosI = griddata2(CenterLon,CenterLat,vlos,LonI,LatI,'nearest');

%% regrid by expanding quads
kount = 0;
for j = 1:nquads
    lons = [min(CornersLon(j,:))+dLon/2.0:dLon:max(CornersLon(j,:))-dLon/2.0]';
    lats = [min(CornersLat(j,:))+dLat/2.0:dLat:max(CornersLat(j,:))-dLat/2.0]';
    for klon = 1:numel(lons)
        lon1 = lons(klon);
        for klat = 1:numel(lats)
            kount = kount+1;
            CenterLon2(kount) = lon1;
            CenterLat2(kount) = lats(klat);
            vlos2(kount) = vlos(j);
        end
    end
end
vlosI = griddata2(CenterLon2,CenterLat2,vlos2,LonI,LatI,'nearest');
figure;
plot(CenterLon2,CenterLat2,'k.');

% write_table(sprintf('%s.lonlatvlos',strrep(qfname,'.out',''))...
%     ,[CenterLon2;CenterLat2;vlos2]...
%     ,'',{'CenterLon2[deg]','CenterLat2[deg]','vlos2'},3,0);

% csvwrite(sprintf('%s_lonlatvlos.csv',strrep(qfname,'.out',''))...
%     ,[colvec(CenterLon2),colvec(CenterLat2),colvec(vlos2)]);

dlmwrite(sprintf('%s_lonlatvlos.csv',strrep(qfname,'.out',''))...
    ,[colvec(CenterLon2),colvec(CenterLat2),colvec(vlos2)]...
    ,'precision','%.8f','delimiter',',');


           

%% project Centers into UTM meters
[CenterXutm,CenterYutm,utmzone] = deg2utm(CenterLat,CenterLon);

%% project corners into UTM meters
CornersXutm = nan(size(CornersLon));
CornersYutm = nan(size(CornersLat));
for i = 1:4
    [CornersXutm(:,i),CornersYutm(:,i),utmzone] = deg2utm(CornersLat(:,i),CornersLon(:,i));
end

%% make figure in geographic coordinates
figure;hold on;axis equal;
patch(CornersLon',CornersLat',vlos*1.e3);
plot(CenterLon,CenterLat,'.k');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title(qfname,'Interpreter','None');
h=colorbar;title(h,'mm/yr');

%% make figure in regridded coordinates
figure;hold on;axis equal;
imagesc(vlosI);

figure;hold on;axis equal;
imagesc([Lonmin:dLon:Lonmax],[Latmin:dLat:Latmax],vlosI);


%% make figure in cartographic UTM coordinates
figure;hold on;axis equal;axis off;
patch(CornersXutm'/1.e3,CornersYutm'/1.e3,vlos*1.e3);
colormap('gray');
% extract pixel values
F=getframe(gcf);

% continue plotting
axis on;
colormap('jet');
plot(CenterXutm/1.e3,CenterYutm/1.e3,'.k');
xlabel('UTM easting [km]');
ylabel('UTM Northing [km]');
title(qfname,'Interpreter','None');
h=colorbar;title(h,'mm/yr');

%% try to read image as an array
figure;
[FX,MAP] = frame2im(F);
FXI=rgb2gray(FX);

%imagesc(F.cdata(:,:,1));
imagesc(FX)
xlabel('UTM easting [km]');
ylabel('UTM Northing [km]');
title(qfname,'Interpreter','None');
h=colorbar;title(h,'mm/yr');






    
