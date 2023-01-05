% test indexing

%A=[11,12;21,22]
A=[1,2,3;-3,-2,-1]

[amin,kmin]=min(colvec(A))

[imin,jmin] = ind2sub(size(A),kmin)