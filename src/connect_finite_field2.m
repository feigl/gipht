function [I,J] = connect_finite_field2(F)
% given a 2-D field with NaN values, return the indices I,J that will
% connect the elements in the field that are finite
% 2018/09/27 Kurt Feigl
[nrows,mcols] = size(F);
npix = nrows*mcols;
% I=zeros(npix^2,1);
% J=zeros(npix^2,1);
%kount = 0;
k1=1:npix;
[i1,j1] = ind2sub([nrows,mcols],k1);  
k2=min([k1+1,npix-1]):npix;
[i2,j2] = ind2sub([nrows,mcols],k2);
K1=find(isfinite(F(i1,j1)));
K2=find(isfinite(F(i2,j2)));
[I1,J1]=ind2sub([nrows,mcols],K1);
[I2,J2]=ind2sub([nrows,mcols],K2);
I=[I1,J1];
J=[I2,J2];

% %                 fprintf(1,'%2d %2d ',i1,j1);
% kount = kount+1;
% I(kount)=i1;
% J(kount)=j1;
% %                 fprintf(1,' %2d %2d\n',i2,j2);
% kount = kount+1;
% I(kount)=i2;
% J(kount)=j2;
% I=I(I>0);
% J=J(J>0);
return
end

