function [I1,J1,I2,J2] = connect_finite_field4(F,nselect)
% given a 2-D field with NaN values, return the indices I,J that will
% connect the elements in the field that are finite
% 2018/09/27 Kurt Feigl
[nrows,mcols] = size(F);
npix = nrows*mcols;
% I=zeros(npix^2,1);
% J=zeros(npix^2,1);
I1=zeros(1,20*nselect);
J1=zeros(1,20*nselect);
I2=zeros(1,20*nselect);
J2=zeros(1,20*nselect);
kount = 0;
%for k1=1:npix
kr1 = randi(npix,[1,20*nselect]);
kr2 = randi(npix,[1,20*nselect]);

for i=1:numel(kr1)
    [i1,j1] = ind2sub([nrows,mcols],kr1(i));
    [i2,j2] = ind2sub([nrows,mcols],kr2(i));
    if isfinite(F(i1,j1)) && isfinite(F(i2,j2)) % && kount < nselect       
        kount = kount+1;
        I1(kount)=i1;
        J1(kount)=j1;     
        I2(kount)=i2;
        J2(kount)=j2;
    end
end
I1=I1(I1>0);
J1=J1(J1>0);
I2=I2(I2>0);
J2=J2(J2>0);
%nselect
I1=I1(1:nselect);
I2=I2(1:nselect);
J1=J1(1:nselect);
J2=J2(1:nselect);

return
end

