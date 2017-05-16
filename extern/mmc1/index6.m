function [m,n]=index6(i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index6.m
% converts from a vector index to a tensor index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if i==1
	m=1;n=1;
elseif i==2
	m=1;n=2;
elseif i==3
	m=1;n=3;
elseif i==4
	m=2;n=2;
elseif i==5
	m=2;n=3;
elseif i==6
	m=3;n=3;
end
