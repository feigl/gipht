function [A_21,A_22]=buildtensors(a)

% builds tensors of up to rank 2 for elementwise multiplication to avoid
% nested for loop evaluations

% Output naming convention is inputvector_## where the first # is the
% tensor order and the second # is the coordinate direction in which the
% elements of the input vector are advanced (i.e. in which the elements are
% unique)
A_11 = a;
A_21=cat(2,A_11',A_11',A_11');
A_22=cat(2,a(1)*ones(3,1),a(2)*ones(3,1),a(3)*ones(3,1));
