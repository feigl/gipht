function [t,x,y,z,u,v,w] = read_abaqus_csv2d(csv_filename)

txzuw = csvread('mogi-test.csv');
t=txzuw(:,1);
x=txzuw(:,2);
z=txzuw(:,3);
u=txzuw(:,4);
v=txzuw(:,5);

