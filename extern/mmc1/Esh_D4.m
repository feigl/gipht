function D4 = Esh_D4(vm,a,x)

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
    I=(acosh(bbar))*4*pi*a(1)*a(2)^2/sqrt(a(1)^2-a(2)^2);
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
    I=(acos(bnonbar))*4*pi*a(1)^2*a(3)/sqrt(a(1)^2-a(3)^2);
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
c2=sum((x.^2)./((a.^2+lambda).^3));
c3=sum((x.^2)./((a.^2+lambda).^4));
    
F=2*x./(a.^2+lambda); 
[F_21,F_22]=buildtensors(F);

if lambda==0
    fderlambda=zeros(1,3);
else
    fderlambda=F/c1;
end
[fderlambda_21,fderlambda_22]=buildtensors(fderlambda);

diagvals=eye(3);
nondiagvals=ones(3)-eye(3);
fderF=nondiagvals.*(1./(a_21.^2+lambda)).*(-F_21.*fderlambda_22)+diagvals.*(1./(a_21.^2+lambda)).*(2-F_21.*fderlambda_22);
fderc1=F./(a.^2+lambda)-2*c2*fderlambda;
[fderc1_21,fderc1_22]=buildtensors(fderc1);
fderc2=F./(a.^2+lambda).^2-3*c3*fderlambda;
[fderc2_21,fderc2_22]=buildtensors(fderc2);

if lambda==0
    sderlambda=zeros(3);
else
    sderlambda=(fderF-fderlambda_21.*fderc1_22)/c1;
end

sderc1=(1./(a_21.^2+lambda)).*(fderF-fderlambda_22.*F_21./(a_21.^2+lambda))-2*(fderc2_22.*fderlambda_21+c2*sderlambda);
fderIfir=ultadelfir_21.*fderlambda_22;

for q=1:3
    for p=1:3
        for r=1:3
            sderF(q,p,r)=-(fderF(q,p)*fderlambda(r)+fderF(q,r)*fderlambda(p)+F(q)*sderlambda(p,r))/(a(q)^2+lambda);
        end
    end
end

zeefir=1./(a.^2+lambda)+0.5*sum(1./(a.^2+lambda));
zeesec=1./(a_21.^2+lambda)+1./(a_22.^2+lambda)+0.5*sum(1./(a.^2+lambda));
for i=1:3
    for j=1:3
        for k=1:3
            sderIfir(i,j,k)=ultadelfir(i)*(sderlambda(j,k)-fderlambda(j)*fderlambda(k)*zeefir(i));
        end
    end
end

for i=1:3
    for j=1:3
        for k=1:3
            fderIsec(i,j,k)=ultadelsec(i,j)*fderlambda(k);
        end
    end
end

for i=1:3;
    for j=1:3
        for k=1:3
            for l=1:3
                sderIsec(i,j,k,l)=ultadelsec(i,j)*(sderlambda(k,l)-fderlambda(k)*fderlambda(l)*zeesec(i,j));
            end
        end
    end
end

for q=1:3
    for p=1:3
        for r=1:3
            if lambda==0
                tderlambda(q,p,r)=0;
            else
                tderlambda(q,p,r)=(-1/c1)*(sderlambda(q,p)*fderc1(r)-sderF(q,p,r)+sderlambda(q,r)*fderc1(p)+fderlambda(q)*sderc1(p,r));
            end
        end
    end
end

%************************************************
%Calculation of V-potentials
%***********************************************

for i=1:3
    for p=1:3
        for q=1:3 
            sderVfir(i,p,q)=-(kdelta(p,q)*Isec(p,i)+x(p)*fderIsec(p,i,q));
        end
    end
end

for i=1:3
    for p=1:3
        for q=1:3
            for r=1:3
                tderVfir(i,p,q,r)=-(kdelta(p,q)*fderIsec(p,i,r)+kdelta(p,r)*fderIsec(p,i,q)+x(p)*sderIsec(p,i,q,r));
            end
        end
    end
end

%*********************************************
%calculation of phi and psi potentials
%*********************************************

%calculation of phi derivatives
for p=1:3
    for q=1:3
        sderphi(p,q)=-(kdelta(p,q)*Ifir(p)+x(p)*fderIfir(p,q));
    end
end

for p=1:3
    for q=1:3
        for r=1:3
            tderphi(p,q,r)=-(kdelta(p,q)*fderIfir(p,r)+kdelta(p,r)*fderIfir(p,q)+x(p)*sderIfir(p,q,r));
        end
    end
end

%*******************
%psi's
%***************

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                foderpsi(i,j,k,l)=kdelta(i,j)*(sderphi(k,l)-a(i)^2*sderVfir(i,k,l))+kdelta(i,k)*(sderphi(j,l)-a(i)^2*sderVfir(i,j,l))+kdelta(i,l)*(sderphi(j,k)-a(i)^2*sderVfir(i,j,k))+x(i)*(tderphi(j,k,l)-a(i)^2*tderVfir(i,j,k,l));
            end
        end
    end
end

%*******************************************
%calculation of D4 
%******************************************
premult1=1/(8*pi*(1-vm));

%calculation of D4 

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                D4(i,j,k,l)=premult1*(foderpsi(k,l,i,j)-2*vm*kdelta(k,l)*sderphi(i,j)-(1-vm)*(sderphi(k,j)*kdelta(i,l)+sderphi(k,i)*kdelta(j,l)+sderphi(l,j)*kdelta(i,k)+sderphi(l,i)*kdelta(j,k)));
            end
        end
    end
end
return