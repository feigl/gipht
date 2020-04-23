function [I,J] = connect_finite_field(F,nselect)
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
for k1=randperm(npix)
    [i1,j1] = ind2sub([nrows,mcols],k1);
    if isfinite(F(i1,j1))
        for k2=min([k1+1,npix-1]):npix
            [i2,j2] = ind2sub([nrows,mcols],k2);
            if isfinite(F(i2,j2)) && kount < nselect
%                 fprintf(1,'%2d %2d ',i1,j1);
                kount = kount+1;
                I(kount)=i1;
                J(kount)=j1;
%                 fprintf(1,' %2d %2d\n',i2,j2);
                kount = kount+1;
                I(kount)=i2;
                J(kount)=j2;
            end
        end
    end
end
I=I(I>0);
J=J(J>0);
return
end

