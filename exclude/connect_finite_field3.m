function [I,J] = connect_finite_field3(F,nselect)
% given a 2-D field with NaN values, return the indices I,J that will
% connect the elements in the field that are finite
% 2018/09/27 Kurt Feigl
[nrows,mcols] = size(F);
npix = nrows*mcols;
% I=zeros(npix^2,1);
% J=zeros(npix^2,1);
I=zeros(nselect,1);
J=zeros(nselect,1);
kount = 0;
%for k1=1:npix
kr1 = randi(npix,[5*nselect,1]);
kr2 = randi(npix,[5*nselect,1]);

for i=1:numel(kr1)
    [i1,j1] = ind2sub([nrows,mcols],kr1(i));
    [i2,j2] = ind2sub([nrows,mcols],kr2(i));
    if isfinite(F(i1,j1)) && isfinite(F(i2,j2)) && kount < nselect
        % fprintf(1,'%2d %2d ',i1,j1);
        kount = kount+1;
        I(kount)=i1;
        J(kount)=j1;
        % fprintf(1,' %2d %2d\n',i2,j2);
        kount = kount+1;
        I(kount)=i2;
        J(kount)=j2;
    end
end
I=I(I>0);
J=J(J>0);
return
end

