%% Comparing GPS measurements to InSAR data
% Elena Reinisch 20180928
% get updated data files with full deployment period
% wget http://geodesy.unr.edu/gps_timeseries/tenv3/NA12/BRD1.NA12.tenv3
%  wget http://geodesy.unr.edu/gps_timeseries/tenv3/NA12/BRDY.NA12.tenv3

% Initialize
clear all; close all;
addpath('Comparison_Tests/GPS')

% load GPS data
BRDY_table = readtable('BRDY.NA12.tenv3.txt', 'ReadVariableNames', 1, 'HeaderLines', 0);
BRD1_table = readtable('BRD1.NA12.tenv3.txt', 'ReadVariableNames', 1, 'HeaderLines', 0);


[x_brdy, y_brdy] = deg2utm(39.781, -119.009);
[x_brd1, y_brd1] = deg2utm(39.809, -119.003);
z_brdy = 1269.092;
z_brd1 = 1247.393 ;


% save data to structure
GPS = struct;
GPS.BRDY = table2struct(BRDY_table, 'toScalar', 1);
GPS.BRD1 = table2struct(BRD1_table, 'toScalar', 1);

% convert dates to datetime format for comparison with InSAR
[tm_y, tm_m, tm_d] = dyear2date(GPS.BRDY.yyyy_yyyy);

GPS.BRDY.datetime = datetime(tm_y, tm_m, tm_d);
[tm_y, tm_m, tm_d] = dyear2date(GPS.BRD1.yyyy_yyyy);
GPS.BRD1.datetime = datetime(tm_y, tm_m, tm_d);

%% re-create plots
% BRD1 - BRDY
[tB1By, tB1, tBY] = intersect(GPS.BRD1.yyyy_yyyy(GPS.BRD1.datetime < '13-Apr-2019'), GPS.BRDY.yyyy_yyyy(GPS.BRDY.datetime < '13-Apr-2019'));
GPS.BRD1_BRDY.tB1BY = tB1By;
GPS.BRD1_BRDY.tB1 = tB1;
GPS.BRD1_BRDY.tBY = tBY;

% calculate relative displacements
dxeast_B1BY = ((GPS.BRD1.x__east_m_(tB1)) - (GPS.BRDY.x__east_m_(tBY)))*1000;
dxnorth_B1BY = ((GPS.BRD1.x_north_m_(tB1)) - (GPS.BRDY.x_north_m_(tBY)))*1000;
dxup_B1BY = ((GPS.BRD1.x____up_m_(tB1)) - (GPS.BRDY.x____up_m_(tBY)))*1000;
sig_east = sqrt([GPS.BRD1.sig_e_m_(tB1).^2 + GPS.BRDY.sig_e_m_(tBY).^2])*1000;
sig_north = sqrt([GPS.BRD1.sig_n_m_(tB1).^2 + GPS.BRDY.sig_n_m_(tBY).^2])*1000;
sig_vert = sqrt([GPS.BRD1.sig_u_m_(tB1).^2 + GPS.BRDY.sig_u_m_(tBY).^2])*1000;

GPS.BRD1_BRDY.dxeast = dxeast_B1BY;
GPS.BRD1_BRDY.dxnorth = dxnorth_B1BY;
GPS.BRD1_BRDY.dxup = dxup_B1BY;
GPS.BRD1_BRDY.sig_east = sig_east;
GPS.BRD1_BRDY.sig_north = sig_north;
GPS.BRD1_BRDY.sig_vert = sig_vert;
GPS.BRD1_BRDY.tdatetime = (GPS.BRD1.datetime(tB1));

% calendar date of starting epoch
[yy mm dd] = dyear2date(tB1By(1))

% fit trend to noisy data
unitv = -1*[0.5559;   -0.1011;    0.8251];
GPS.BRDY.range = [GPS.BRDY.x__east_m_(tBY), GPS.BRDY.x_north_m_(tBY), GPS.BRDY.x____up_m_(tBY)]*unitv;
GPS.BRD1.range = [GPS.BRD1.x__east_m_(tB1), GPS.BRD1.x_north_m_(tB1), GPS.BRD1.x____up_m_(tB1)]*unitv;
%%
[P,S] = polyfit(GPS.BRD1.yyyy_yyyy(tB1),GPS.BRD1.range,3)
figure; plot(GPS.BRD1.yyyy_yyyy(tB1), GPS.BRD1.range, 'ro-')
hold on
%tplot = linspace(GPS.BRD1.yyyy_yyyy(tB1(1)), GPS.BRD1.yyyy_yyyy(tB1(end)), 100)
brd1_poly = polyval(P, GPS.BRD1.yyyy_yyyy(tB1));
plot(GPS.BRD1.yyyy_yyyy(tB1), brd1_poly, 'b-')

GPS.BRD1.range_poly = brd1_poly;

[P,S] = polyfit(GPS.BRDY.yyyy_yyyy(tBY),GPS.BRDY.range,3)
figure; plot(GPS.BRDY.yyyy_yyyy(tBY), GPS.BRDY.range, 'ro-')
hold on
brdy_poly = polyval(P, GPS.BRDY.yyyy_yyyy(tBY));
plot(GPS.BRDY.yyyy_yyyy(tBY), brdy_poly, 'b-')

GPS.BRDY.range_poly = brdy_poly;


save('BRD1_BRDY_new_full2.mat', 'GPS')
return
% plot
figure;
errorbar((tB1By), dxeast_B1BY-dxeast_B1BY(1), sig_east)
% hold on;
% plot((tB1By), dxeast_B1BY, 'ko')
title('BRD1 - BRDY Eastward Displacement Relative to March 5, 2016')
xlabel('year')
ylabel('mm')
axis([2016, 2017, -9, 10])
printeps_e('BRD1_BRDY_east.eps')

figure;
errorbar((tB1By), dxnorth_B1BY-dxnorth_B1BY(1), sig_north)
% hold on;
% plot((tB1By), dxnorth_B1BY, 'ko')
title('BRD1 - BRDY Northward Displacement Relative to March 5, 2016')
xlabel('year')
ylabel('mm')
axis([2016, 2017, -9, 13])
printeps_e('BRD1_BRDY_north.eps')


figure;
errorbar((tB1By), dxup_B1BY-dxup_B1BY(1), sig_vert)
% hold on;
% plot((tB1By), dxup_B1BY, 'ko')
title('BRD1 - BRDY Vertical Displacement Relative to March 5, 2016')
xlabel('year')
ylabel('mm')
axis([2016, 2017, -15, 25])
printeps_e('BRD1_BRDY_vert.eps')

