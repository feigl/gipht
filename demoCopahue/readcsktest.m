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


close all;
clear all;

nbytes = 59885568
mcol = 2432
nrow = 3078

lat1 = -37.562215048
dlat = -0.000277777
lat2 = lat1 + nrow * dlat

lon1 = -71.528334654;
dlon = 0.000277777
lon2 = lon1 + mcol * dlon;


nbytes/4/2
nrow*mcol

bil_file_name = '20111023_20111022/filt_topophase.unw.geo'
[amp,pha] = read_bil(bil_file_name,nrow,mcol);

% ibad = find(amp>1.e7);
% amp(ibad)=nan;

amp = log(amp);

figure; hist(colvec(amp),255);
max(max(double(amp)))

figure;hold on;
imagesc(lon1:dlon:lon2,lat1:dlat:lat2,double(amp));
colormap('gray');colorbar
title(sprintf('amplitude %s',bil_file_name),'Interpreter','none');

figure;hold on;axis equal;axis ij;
imagesc(lon1:dlon:lon2,lat1:dlat:lat2,double(pha));
colormap('jet');colorbar
title(sprintf('phase %s',bil_file_name),'Interpreter','none');
