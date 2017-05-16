function [S4,PIvector]=Eshint(vm,a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eshint.m
% Calculates the internal Eshelby tensors (see
% Eshelby '57 eqn. 3.8) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% Changes made
% error checking
% change all var names to get rid of or1in1
% make the case decisions dependent on 1e-6 * dimension, instead of the absolute size 1e-6
% fix the case decision logic
% fix case and statements to use logical and (&&) rather than elementwise and (&)
% general formatting
% verified and corrected I formulations from Mura

%% Changes to do
% Finish verification of functions starting in line 97

%******************************************************************%
%Calculation of I's
%******************************************************************%

if a(1)<=0 || a(2)<=0 || a(3)<=0
error('Ellipsoid dimensions (a) must be positive')
end

if abs(a(1)-a(2))<(1e-6*a(1)) && abs(a(2)-a(3))<(1e-6*a(1))  % checks that geometric mean of ellipsoid dimensions is not more than 1e-6 different from first dimension
    % Spherical Case
    Ifir=(4/3)*pi*ones(1,3);
    Isec=(4/5)*pi*a(1)^2*ones(3);
    
elseif (a(1)-a(2))>(1e-6*a(1)) && abs(a(3)-a(2))<(1e-6*a(1))
    % Prolate Spheriod Case
    rat=a(1)/a(3);	
    
    Ifir=zeros(1,3);
    Ifir(2)=(2*pi*a(1)*a(3)^2/((a(1)^2-a(3)^2)^(3/2)))*(rat*(rat^2-1)^(1/2)-acosh(rat));
    Ifir(3)=Ifir(2);
    Ifir(1)=4*pi-2*Ifir(2);
    
    Isec=zeros(3);
    Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)^2-a(2)^2);
    Isec(1,3)=Isec(1,2);
    Isec(2,1)=Isec(1,2);
    Isec(3,1)=Isec(1,3);
    Isec(1,1)=(4*pi/a(1)^2-2*Isec(1,2))/3;
    Isec(2,3)=pi/(a(2)^2)-(Ifir(2)-Ifir(1))/(4*(a(1)^2-a(2)^2));
    Isec(3,2)=Isec(2,3);
    Isec(2,2)=Isec(2,3);
    Isec(3,3)=Isec(2,3);

elseif abs(a(1)-a(2))<(1e-6*a(1)) && (a(2)-a(3))>(1e-6*a(2))
    % Oblate Spheriod Case
    rat=a(3)/a(1);	
    
    Ifir=zeros(1,3);
    Ifir(1)=(2*pi*a(1)^2*a(3)/((a(1)^2-a(3)^2)^(3/2)))*(acos(rat)-rat*(1-rat^2)^(1/2));
    Ifir(2)=Ifir(1);
    Ifir(3)=4*pi-2*Ifir(1);
    
    Isec=zeros(3);
    Isec(1,3)=(Ifir(1)-Ifir(3))/(a(3)^2-a(1)^2);
    Isec(3,1)=Isec(1,3);
    Isec(2,3)=Isec(1,3);
    Isec(3,2)=Isec(2,3);
    Isec(1,2)=pi/a(1)^2-Isec(1,3)/4;
    Isec(2,1)=Isec(1,2);
    Isec(1,1)=Isec(1,2);
    Isec(2,2)=Isec(1,2);
    Isec(3,3)=(4*pi/a(3)^2-2*Isec(1,3))/3;
    
else
    % Triaxial Ellipsoid Case    
    theta=asin((1-(a(3)/a(1))^2)^(1/2)); % amplitude
    % k=((a(1)^2-a(2)^2)/(a(1)^2-a(3)^2))^(1/2); % the elliptic modulus
    m=(a(1)^2-a(2)^2)/(a(1)^2-a(3)^2); % m=k^2 is the parameter;
    [F,E,Z]=elliptic12(theta,m); %this sets the tolerance to eps, add a third argument to set to a larger tol
    % Mura 11.17
    Ifir=zeros(1,3);
    Ifir(1)=(4*pi*prod(a)/((a(1)^2-a(2)^2)*sqrt(a(1)^2-a(3)^2)))*(F-E);
    Ifir(3)=(4*pi*prod(a)/((a(2)^2-a(3)^2)*sqrt((a(1)^2-a(3)^2))))*(a(2)*sqrt((a(1)^2-a(3)^2))/(a(1)*a(3))-E);
    Ifir(2)=4*pi-Ifir(1)-Ifir(3);
    
    Isec=zeros(3);
    Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)^2-a(2)^2);
    Isec(2,3)=(Ifir(3)-Ifir(2))/(a(2)^2-a(3)^2);
    Isec(3,1)=(Ifir(1)-Ifir(3))/(a(3)^2-a(1)^2);
    Isec(2,1)=Isec(1,2);
    Isec(3,2)=Isec(2,3);
    Isec(1,3)=Isec(3,1);
    Isec(1,1)=(4*pi/a(1)^2-Isec(1,2)-Isec(1,3))/3;
    Isec(2,2)=(4*pi/a(2)^2-Isec(2,3)-Isec(2,1))/3;
    Isec(3,3)=(4*pi/a(3)^2-Isec(3,1)-Isec(3,2))/3;    
end

denom=8*pi*(1-vm);

S1111=(3*a(1)^2*Isec(1,1)+(1-2*vm)*Ifir(1))/denom;
S2222=(3*a(2)^2*Isec(2,2)+(1-2*vm)*Ifir(2))/denom;
S3333=(3*a(3)^2*Isec(3,3)+(1-2*vm)*Ifir(3))/denom;

S1122=(a(2)^2*Isec(1,2)-(1-2*vm)*Ifir(1))/denom;
S2233=(a(3)^2*Isec(2,3)-(1-2*vm)*Ifir(2))/denom;
S3311=(a(1)^2*Isec(3,1)-(1-2*vm)*Ifir(3))/denom;

S1133=(a(3)^2*Isec(1,3)-(1-2*vm)*Ifir(1))/denom;
S2211=(a(1)^2*Isec(2,1)-(1-2*vm)*Ifir(2))/denom;
S3322=(a(2)^2*Isec(3,2)-(1-2*vm)*Ifir(3))/denom;

S1212=((a(1)^2+a(2)^2)*Isec(1,2)+(1-2*vm)*(Ifir(1)+Ifir(2)))/(2*denom);
S2323=((a(2)^2+a(3)^2)*Isec(2,3)+(1-2*vm)*(Ifir(2)+Ifir(3)))/(2*denom);
S3131=((a(3)^2+a(1)^2)*Isec(3,1)+(1-2*vm)*(Ifir(3)+Ifir(1)))/(2*denom);
S1313=S3131;

S4=[S1111,  0,      0,      S1122,  0,      S1133;
    0,      2*S1212,  0,      0,      0,      0;
    0,      0,      2*S1313,  0,      0,      0;
    S2211,  0,      0,      S2222,  0,      S2233;
    0,      0,      0,      0,      2*S2323,  0;
    S3311,  0,      0,      S3322,  0,      S3333];

PI3131=(Ifir(1)-Ifir(3))/(8*pi);
 
PI1212=(Ifir(2)-Ifir(1))/(8*pi);

PI2323=(Ifir(3)-Ifir(2))/(8*pi);
 
PI1313=-PI3131;
PI2121=-PI1212;
PI3232=-PI2323;

PIvector=[2*PI3232;2*PI1313;2*PI2121];
