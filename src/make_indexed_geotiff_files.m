close all
clear all

addpath(genpath('/Users/feigl/gipht/src'),'-begin');
addpath(genpath('/Users/feigl/gipht/extern'),'-begin');
addpath(genpath('/Users/feigl/gipht/utils'),'-begin');

% on porotomo
% cd /t31/insar/TSX/T53/brady/Best_Pairs/intf/In20111224_20121027
% gmt grdhisteq display_amp_utm.grd -Gdisplay_amp_eq_utm.grd -V
% source ~ebaluyut/setup.sh
% which grd2geotiff.csh
% view /t31/ebaluyut/scratch/TEST_GMTSAR/gmtsar/bin/grd2geotiff.csh
% grd2geotiff.csh display_amp_eq_utm.grd

% make a color table with 2^16 values
map16 = nan(2^16,3);
map8  = hsv;
k=0;
for i=1:64
    for j=1:256*4
        k = k+1;
        map16(k,1) = map8(i,1);
        map16(k,2) = map8(i,2);
        map16(k,3) = map8(i,3);
    end
end
% make first entry gray
map16(1:64,1) = 0.5;
map16(1:64,2) = 0.5;
map16(1:64,3) = 0.5;


%infilenames = {'display_amp_utm.tif','drhomaskd_utm.tif','phasefilt_utm.tif','display_amp_eq_utm.tif'}
infilenames = {'drhomaskd_utm.tif','phasefilt_utm.tif','phasefilt_mask_ll.tif'}

for i=1:numel(infilenames)
    infilename1 = infilenames{i};
    outfilename = strrep(infilename1,'.tif','_ind.tif')
    lut = geotiff_float_to_16bit(infilename1,outfilename,map16);
end

% figure
% Amp = imread('display_amp_utm_ind.tif');
% imshow(Amp,map)
% colormap(gray)
% 
% 
% figure
% Aeq = imread('display_amp_eq_utm_ind.tif');
% imshow(Aeq,map)
% colormap(gray)

%% read phase file
figure;hold on;
[Pha,map] = imread('phasefilt_utm_ind.tif');
imshow(Pha,map); hold on;
[inan2,jnan2] =find(Pha==1);
%plot(jnan2,inan2,'k.');
colorbar



