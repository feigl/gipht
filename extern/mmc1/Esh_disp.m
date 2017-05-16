function u=Esh_disp(vm,a,x,eigen)

% todo search for todos in function
% are the case statements supposed to be exact or with a tolerance like in eshint
% get rid of all vars not used

%******************************************************************%
%Calculation of F and E integrals
%******************************************************************%

% this subroutines finds the largest positive root of
% x(1)^2/(a(1)+lambda) + x(2)^2/(a(2)+lambda) + x(3)^2/(a(3)+lambda) = 1
% (Mura 11.37) for the exterior point x and elliopsoid dimensions a.  When 
% expanded and like terms in lambda are collected, the coefficients of 
% lambda^3, ^2, etc. are as below

coef3=1; % coefficient of lambds^3 term
coef2=a(1)^2+a(2)^2+a(3)^2-(x(1)^2+x(2)^2+x(3)^2); % coefficient of lambds^2 term
coef1=a(1)^2*a(2)^2+a(1)^2*a(3)^2+a(2)^2*a(3)^2-((a(2)^2+a(3)^2)*x(1)^2+(a(1)^2+a(3)^2)*x(2)^2+(a(1)^2+a(2)^2)*x(3)^2); % coefficient of lambds term
coef0=a(1)^2*a(2)^2*a(3)^2-(a(2)^2*a(3)^2*x(1)^2+a(1)^2*a(3)^2*x(2)^2+a(1)^2*a(2)^2*x(3)^2); % coefficient of constant term
poly=[coef3,coef2,coef1,coef0]; % matlab polynomial format
lambda=0; % initialize lambda to zero

if x(1)^2/a(1)^2+x(2)^2/a(2)^2+x(3)^2/a(3)^2>1 % if x is exterior point set
    % lambda to the largest positive real root, otherwise lambda=0
    lambdaroots=roots(poly); % store the roots of the cubic equation
    for i=1:3 % find the largest positive real root
        if isreal(lambdaroots(i)) && lambdaroots(i)>lambda
            lambda=lambdaroots(i);
        end
    end
end

theta=asin(((a(1)^2-a(3)^2)/(a(1)^2+lambda))^(1/2)); % the amplitude
% todo this argument was taken from the previous code (with the lambda) and
% modified with the arcsin.  need to see if can get here via Gradshteyn and
% Ryzhik from Mura 11.36
% k=((a(1)^2-a(2)^2)/(a(1)^2-a(3)^2))^(1/2); % the elliptic modulus
m=(a(1)^2-a(2)^2)/(a(1)^2-a(3)^2); % m=k^2 is the parameter
[F,E,Z]=elliptic12(theta,m); %this sets the tolerance to eps, add a third argument to set to a larger tol

%******************************************************************%
%Calculation of I's
%******************************************************************%

if a(1)==a(2) && a(1)==a(3)
    % Spherical Case
    del=sqrt(prod(a.^2+lambda));
    % can simplify to del3=sqrt((a(1)^2+lambda)^3) for sphere
    Ifir=(4/3)*pi*a(1)^3/(a(1)^2+lambda)^(3/2)*ones(1,3);
    Isec=(4/5)*pi*a(1)^3/(a(1)^2+lambda)^(1/2)*ones(3);  % todo: i changed the 5/2 to 1/2 to make units right--not sure if correct

elseif a(1)>a(2)&& a(3)==a(2)
    %fprintf('Prolate case..\n')
    
    del=sqrt((a(1)^2+lambda)*(a(2)^2+lambda)*(a(3)^2+lambda));
    bbar=sqrt(a(1)^2+lambda)/sqrt(a(3)^2+lambda);
    dbar=sqrt(a(1)^2-a(3)^2)/sqrt(a(3)^2+lambda);
    %I=(acosh(bbar))*4*pi*a(1)*a(2)^2/sqrt(a(1)^2-a(2)^2);
    Ifir(1)=4*pi*a(1)*a(2)^2*(acosh(bbar)-dbar/bbar)/(a(1)^2-a(2)^2)^1.5;
    Ifir(2)=2*pi*a(1)*a(2)^2*(-acosh(bbar)+dbar*bbar)/(a(1)^2-a(2)^2)^1.5;
    Ifir(3)=Ifir(2);

    Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)^2-a(2)^2);
    Isec(1,3)=Isec(1,2);
    Isec(2,1)=Isec(1,2);
    Isec(3,1)=Isec(1,3);
    Isec(2,3)=pi*prod(a)/((a(3)^2+lambda)*del)-Isec(1,3)/4;
    Isec(3,2)=Isec(2,3);
    Isec(1,1)=((4*pi*prod(a))/((a(1)^2+lambda)*del)-Isec(1,2)-Isec(1,3))/3;
    Isec(2,2)=Isec(2,3);
    Isec(3,3)=Isec(2,3);
    
elseif a(1)==a(2)&& a(2)>a(3)
    %fprntf('Oblate case...\n')
    del=sqrt((a(1)^2+lambda)*(a(2)^2+lambda)*(a(3)^2+lambda));
    bnonbar=sqrt(a(3)^2+lambda)/sqrt(a(1)^2+lambda);
    dnonbar=sqrt(a(1)^2-a(3)^2)/sqrt(a(1)^2+lambda);
    %I=(acos(bnonbar))*4*pi*a(1)^2*a(3)/sqrt(a(1)^2-a(3)^2);
    Ifir(1)=2*pi*a(1)^2*a(3)*(acos(bnonbar)-dnonbar*bnonbar)/(a(1)^2-a(3)^2)^1.5;
    Ifir(2)=Ifir(1);
    Ifir(3)=4*pi*prod(a)/del-2*Ifir(1);
    
    
    Isec(1,3)=(Ifir(3)-Ifir(1))/(a(1)^2-a(3)^2);
    Isec(3,1)=Isec(1,3);
    Isec(2,3)=Isec(1,3);
    Isec(3,2)=Isec(2,3);
    
    Isec(1,1)=pi*prod(a)/((a(1)^2+lambda)*del)-Isec(1,3)/4;
    Isec(1,2)=Isec(1,1);
    Isec(2,1)=Isec(1,2);
    
    Isec(2,2)=Isec(1,1);
    Isec(3,3)=((4*pi*prod(a))/((a(3)^2+lambda)*del)-Isec(1,3)-Isec(2,3))/3;
else
    %fprintf('triaxial ellipsoid case ..\n')
    del=sqrt((a(1)^2+lambda)*(a(2)^2+lambda)*(a(3)^2+lambda));
    I=4*pi*prod(a)*F/sqrt(a(1)^2-a(3)^2);
    Ifir(1)=I*(1-E/F)/(a(1)^2-a(2)^2); 
    Ifir(2)=4*pi*prod(a)*(E*sqrt(a(1)^2-a(3)^2)/((a(1)^2-a(2)^2)*(a(2)^2-a(3)^2))-F/((a(1)^2-a(2)^2)*sqrt(a(1)^2-a(3)^2))-(1/(a(2)^2-a(3)^2))*sqrt((a(3)^2+lambda)/((a(1)^2+lambda)*(a(2)^2+lambda))));
    Ifir(3)=4*pi*prod(a)/del-Ifir(1)-Ifir(2);
    Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)^2-a(2)^2);
    Isec(2,1)=Isec(1,2);
    Isec(1,3)=(Ifir(3)-Ifir(1))/(a(1)^2-a(3)^2);
    Isec(3,1)=Isec(1,3);
    Isec(2,3)=(Ifir(3)-Ifir(2))/(a(2)^2-a(3)^2);
    Isec(3,2)=Isec(2,3);
    Isec(1,1)=((4*pi*prod(a))/((a(1)^2+lambda)*del)-Isec(1,2)-Isec(1,3))/3;
    Isec(2,2)=((4*pi*prod(a))/((a(2)^2+lambda)*del)-Isec(1,2)-Isec(2,3))/3;
    Isec(3,3)=((4*pi*prod(a))/((a(3)^2+lambda)*del)-Isec(1,3)-Isec(2,3))/3;
end

%*************************************************************************************************
%I derivatives
%*************************************************************************************************

[a_21,a_22]=buildtensors(a);
ultadelfir=-2*pi*prod(a)./((a.^2+lambda)*del);
[ultadelfir_21,ultadelfir_22]=buildtensors(ultadelfir);
ultadelsec=-2*pi*prod(a)./((a_21.^2+lambda).*(a_22.^2+lambda)*del);

% derivatives of lambda
c1=sum((x.^2)./((a.^2+lambda).^2));
F=2*x./(a.^2+lambda); 

if lambda==0
    fderlambda=zeros(1,3);
else
    fderlambda=F/c1;
end
[fderlambda_21,fderlambda_22]=buildtensors(fderlambda);

fderIfir=ultadelfir_21.*fderlambda_22;
for i=1:3
    for j=1:3
       for k=1:3
           fderIsec(i,j,k)=ultadelsec(i,j)*fderlambda(k);
       end
   end
end

% calculate phi derivative
fderphi=-x.*Ifir;

% psi's
for i=1:3
    for j=1:3
        for l=1:3            
                tderpsi(i,j,l)=-kdelta(i,j)*x(l)*(Ifir(l)-a(i)^2*Isec(i,l))-x(i)*x(j)*(fderIfir(j,l)-a(i)^2*fderIsec(i,j,l))-(kdelta(i,l)*x(j)+kdelta(j,l)*x(i))*(Ifir(j)-a(i)^2*Isec(i,j));
        end
    end
end

premult1=1/(8*pi*(1-vm));

%calculate disp
eigenM = [eigen(1:3)';eigen([2 4 5])'; eigen([3 5 6])'];
diag = eigenM(1,1)+eigenM(2,2)+eigenM(3,3); 
u = (premult1*(tprod(tderpsi,[1 -1 -2],eigenM,[-1 -2])-2*vm*diag*fderphi'-4*(1-vm)*(eigenM*fderphi')))';
return
end