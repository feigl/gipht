
% https://www.mathworks.com/matlabcentral/answers/7600-saving-multiple-vectors-with-different-lengths-in-one-matrix
clear all

v1 = [1 2 3];
v2 = [4 5];
% Make padded array.  Could use rows or columns...
% Or use a cell array.
Trees{1} = v1;
Trees{2} = v2;
for i=1:2
    fprintf(1,'%2d ',Trees{i});
    fprintf(1,'\n');
end