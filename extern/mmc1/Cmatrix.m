function matr=Cmatrix(Cm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cmatrix.m
%this function converts the 4th order isotropic stiffness tensor into 6 x 6 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

clear i;
clear j;
clear m;
clear n;
clear p;
clear q;

for i=1:6
for j=1:6
   
    
   [m,n]=index6(i);
   [p,q]=index6(j);

if j==2
matr(i,j)=Cm(m,n,p,q)+Cm(m,n,q,p);
elseif j==3
matr(i,j)=Cm(m,n,p,q)+Cm(m,n,q,p);
elseif j==5
matr(i,j)=Cm(m,n,p,q)+Cm(m,n,q,p);
else
matr(i,j)=Cm(m,n,p,q);
end
end
end

clear i;
clear j;
clear m;
clear n;
clear p;
clear q;