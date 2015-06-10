function [t,x,z,u,w] = read_abaqus_csv2d(csv_filename)
%function [t,x,z,u,w] = read_abaqus_csv2d(csv_filename)
% read a Comma Separated Values files, as output by Abaqus and translated
% by convert_abaqus_dat_to_csv
% 
% 
% Example:
% 
%     [t,x,z,u,w] = read_abaqus_csv2d('mogi-test.csv');
%     quiver(x,z,u,w);
% 
% 
% 20130522 Kurt Feigl

txzuw = csvread(csv_filename);
t=txzuw(:,1);
x=txzuw(:,2);
z=txzuw(:,3);
u=txzuw(:,4);
w=txzuw(:,5);


return

