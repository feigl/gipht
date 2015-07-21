% hengill.geology.wisc.edu% pwd
% /data/Copahue/CSK/20111023_20111022
% hengill.geology.wisc.edu% ls -ltr *
% -rw-r--r--. 1 feigl feigl     2504 Dec  3  2013 filt_topophase.unw.geo.xml
% -rw-r--r--. 1 feigl feigl 59885568 Dec  3  2013 filt_topophase.unw.geo

%         <property name="startingValue">
%             <value>-71.528334654</value>
%             <doc>{'doc': 'Starting value of the coordinate.'}</doc>
%             <units>{'units': 'degree'}</units>
%         </property>
%         <property name="delta">
%             <value>0.000277777</value>
%             <doc>{'doc': 'Coordinate quantization.'}</doc>
%             <units>{}</units>
%         </property>
%         <property name="size">
%             <value>3078</value>
%             <doc>{'doc': 'Coordinate size.'}</doc>
%         </property>
%     </component>
%     <property name="NUMBER_BANDS">
%         <value>2</value>
%     </property>
%     <component name="Coordinate2">
%         <factorymodule>isceobj.Image</factorymodule>
%         <factoryname>createCoordinate</factoryname>
%         <doc>Second coordinate of a 2D image (length).</doc>
%         <property name="startingValue">
%             <value>-37.562215048</value>
%             <doc>{'doc': 'Starting value of the coordinate.'}</doc>
%             <units>{'units': 'degree'}</units>
%         </property>
%         <property name="delta">
%             <value>-0.000277777</value>
%             <doc>{'doc': 'Coordinate quantization.'}</doc>
%             <units>{}</units>
%         </property>
%         <property name="size">
%             <value>2432</value>
%             <doc>{'doc': 'Coordinate size.'}</doc>

% filt_topophase.unw.geo.xml:        <value>l</value>
% filt_topophase.unw.geo.xml:        <value>write</value>
% filt_topophase.unw.geo.xml:        <value>DEM-flattened interferogram orthorectified to an equi-angular latitude, longitude grid</value>
% filt_topophase.unw.geo.xml:        <value>FLOAT</value>
% filt_topophase.unw.geo.xml:        <value>unw</value>
% filt_topophase.unw.geo.xml:        <value>filt_topophase.unw.geo</value>
% filt_topophase.unw.geo.xml:            <value>-71.528334654</value>
% filt_topophase.unw.geo.xml:            <doc>{'doc': 'Starting value of the coordinate.'}</doc>
% filt_topophase.unw.geo.xml:            <value>0.000277777</value>
% filt_topophase.unw.geo.xml:            <value>3078</value>
% filt_topophase.unw.geo.xml:        <value>2</value>
% filt_topophase.unw.geo.xml:            <value>-37.562215048</value>
% filt_topophase.unw.geo.xml:            <doc>{'doc': 'Starting value of the coordinate.'}</doc>
% filt_topophase.unw.geo.xml:            <value>-0.000277777</value>
% filt_topophase.unw.geo.xml:            <value>2432</value>
% filt_topophase.unw.geo.xml:        <value>3078</value>
% filt_topophase.unw.geo.xml:        <value>2432</value>
% filt_topophase.unw.geo.xml:        <value>BIL</value>
% filt_topophase.unw.geo.xml:        <value>Release: 1.5.01, svn-1191, 20131028. Current: svn-1348M.</value>



close all;
clear all;
format long

nbytes = 59885568
% mcol = 2432
% nrow = 3078

nrow = 2432
mcol = 3078

lat1 = -37.562215048
dlat = -0.000277777
lat2 = lat1 + (nrow-1) * dlat

lon1 = -71.528334654;
dlon = 0.000277777
lon2 = lon1 + (mcol-1) * dlon;


nbytes/4/2
nrow*mcol

bil_file_name = '20141217_20130402/filt_topophase.unw.geo'
[ampB,phaB] = read_bil(bil_file_name,nrow,mcol);

amp_file_name = '20141217_20130402/ampli.grd'
[lonA,latA,ampG] = grdread2(amp_file_name);
if abs(min(lonA) - min([lon1,lon2])) > abs(dlon)
    min(lonA)
    min([lon1,lon2])
    abs(dlon)
   error('1')
end
if abs (min(latA) - min([lat1,lat2])) > abs(dlat)
    min(latA)
    min([lat1,lat2])
    abs(min(latA) - min([lat1,lat2]))
    abs(min(latA) - min([lat1,lat2]))/dlat
    error('2')
end

pha_file_name = '20141217_20130402/phase.grd'
[lonP,latP,phaG] = grdread2(pha_file_name);
if abs (min(lonP) - min([lon1,lon2])) > abs(dlon)
    min(lonP) 
    min([lon1,lon2])
    abs(dlon)
   error('3')
end
if abs (min(latP) - min([lat1,lat2])) > abs(dlat)
    min(latP)
    min([lat1,lat2])
    abs(dlat)
    error('4')
end

coh_file_name = '20141217_20130402/coher.grd'
[lonC,latC,cohG] = grdread2(coh_file_name);
if abs (min(lonC) - min([lon1,lon2])) > abs(dlon)
    min(lonC) 
    min([lon1,lon2])
    abs(dlon)
   error('3')
end
if abs (min(latC) - min([lat1,lat2])) > abs(dlat)
    min(latC)
    min([lat1,lat2])
    abs(dlat)
    error('4')
end


% ibad = find(ampB>1.e7);
% ampB(ibad)=nan;

ampB = log(double(ampB));
ampG = log(double(ampG));

% transfer NaNs from amplitude to phase
ibad = find(isfinite(ampG)==0);
phaG(ibad) = nan;

ilow = find(cohG < 0.7);
phaG(ilow) = nan;

% wrap on 31 mm wavelength
phaG = rwrapm(2*pi*(-1.0*phaG/31/2. - 31/4.));

phabyte = 256.*phaG/2./pi;
phabyte_file_name = strrep(pha_file_name,'.grd','.pha');
%write_pha(phabyte_file_name,phabyte);
write_pha(phabyte_file_name,flipud(phabyte));
phabyte2 = double(read_pha(phabyte_file_name,mcol));
izero = find(abs(phabyte2) < 1);
phabyte2(izero) = nan;


% figure; hist(colvec(ampB),255);
% max(max(double(ampB)))

figure;
hist(colvec(double(phabyte2)),256);
title('phabyte2');

figure;
hist(colvec(double(cohG(cohG>0.01))),256);
title('cohG');

figure;hold on;
imagesc(lonC,latC,double(cohG));
colormap('jet');cmapblacknan;colorbar;
title(sprintf('coherence %s',coh_file_name),'Interpreter','none');

figure;hold on;
imagesc(lon1:dlon:lon2,lat1:dlat:lat2,double(ampB));
colormap('gray');colorbar;
title(sprintf('amplitude %s',bil_file_name),'Interpreter','none');

figure;hold on;
imagesc(lon1:dlon:lon2,lat1:dlat:lat2,double(phaB));
colormap('jet');cmapblacknan;colorbar;
title(sprintf('phase %s',bil_file_name),'Interpreter','none');

figure;hold on;
imagesc(lonA,latA,double(ampG));
colormap('gray');colorbar;
title(sprintf('amplitude %s',amp_file_name),'Interpreter','none');

figure;hold on;
imagesc(lonP,latP,double(phaG));
colormap('jet');cmapblacknan;colorbar;
title(sprintf('phase %s',pha_file_name),'Interpreter','none');

figure;hold on;
imagesc(lonP,latP,double(phabyte2));
colormap('jet');cmapblacknan;colorbar;
title(sprintf('phase %s',phabyte_file_name),'Interpreter','none');




