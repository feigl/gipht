function [mest,psig,pest,mse,res] = geostatistical_inversion1(A,Q,b,AR,Xb,zcol)
%% perform geostatistical inversion
% inputs:
%     A == design matrix
%     R == data covariance matrix
%     Q == prior model covariance matrix
%     b == data
%     Xb = X matrix for beta 
%     zcol related to trend

% outputs:
%     mest         == Mode of estimated model parameter values
%     psig         == diagonal elments of posterior model covariance
%     pest         == posterior?
%     mse          == mean squared error (square of residuals)
    
% References

% Reinisch, E. C., M. Cardiff, and K. L. Feigl (2018), Characterizing
% Volumetric Strain at Brady Hot Springs, Nevada, USA Using Geodetic Data,
% Numerical Models, and Prior Information, Geophysical Journal
% International, 1501?1513. http://dx.doi.org/10.1093/gji/ggy347 

% Kitanidis, P. K. (2007), On stochastic inverse modeling, Geophysical
% Monograph-American Geophysical Union, 171, 19.

% Cardiff, M. (2020), Bayesian Approaches to Inverse Problems, edited.

%% verify sizes
[ndata,mparams] = size(A)


[nAR,mAR] = size(AR)
% if nAR ~= ndata || mAR ~= ndata
%     nAR
%     mAR
%     ndata
%     error(sprintf('Data covariance matrix R must be %d by %d\n',ndata,ndata));
% end

[nQ,mQ] = size(Q)
if nQ ~= mparams || mQ ~= mparams
    nQ
    mQ
    mparams
    error(sprintf('Model covariance matrix Q must be %d by %d\n',mparams,mparams));
end

if nAR > 1 && mAR > 1
    include_AR = 1;
else
    include_AR = 0;
end

disp('Start defining Ggeo and Dgeo')

if include_AR == 0
    % for no data covariance
    % multiply by A and Xb to decrease size
    Ggeo = [(A'*A*Q)*A'+Xb*(A*Xb)', A'*(A*Xb)+Xb*(zeros(size(zcol)))];
else
    % for when we have data covariance
    % multiply by A and Xb to decrease size
    C = zeros(zcol);
    C(:,2) = 1; % constrain zero mean to non def region
    Ggeo = [(A'*A*Q)*A'+AR+Xb*(A*Xb)', A'*(A*Xb)+Xb*C];
    %clear AR
end
%clear Q
disp('Ggeo defined')

%load('b.mat')
Dgeo = A'*b+Xb*(zeros(zcol, 1));
%clear b
disp('Dgeo defined')

%clearvars -except A AR ARfac Ggeo Dgeo include_AR V0 vfactor nf depthsz centx centy dx dy datx daty nrowsA ncolsA datx_vec daty_vec grdx grdy dz dt sig2 mu bi ncells lx ly irefcount zcol topz


 % from regular ls with full dVdt3=unit(-23107,'m^3/year'), from 4Okada -3.1624e+04
 % from regular ls with < 500, DV =  -37107.3 m^3
% from Ali et al. (2016) result for In20130513_20140511: 2013.3616 2014.3562 -19244.912 3863.12 TSX T53_32785_38296 0.898067 0.906401 0.9946
%pest = pinv(Ggeo'*Ggeo)*Ggeo'*Dgeo; %doulbe version: DV = -376305 m^3 %pinv(Ggeo'*Ggeo,  1e-18)*Ggeo'*Dgeo; % DV =  -37236.9 m^3
%load('Q.mat')
disp('Start inversion')
pest = pinv((Ggeo))*(Dgeo);

%clear Dgeo

disp('inversion completed')
% load('A.mat')
% load('Xb.mat')
% load('Q.mat')

disp('finding posterior model covariance')

% store only diagonal elements
psig = diag(Q - [A*Q; Xb']'*pinv(Ggeo)*[A', Xb]*[A*Q; Xb']);

disp('defining beta and xi')
%clear Ggeo
%load('b.mat')
% beta = pest(end)
% xi = pest(1:end-1);
beta = pest(end-1:end);
xi = pest(1:end-2);

% mode 
disp('Mode of estimated parameter values')
mest = Xb*beta + Q*A'*xi;

% compare to true model parameters
if include_AR == 0
    mse = (b-A*mest)'*(b-A*mest)
else
    %mse0 = (b-A*mest)'*(b-A*mest)
    res = (b-A*mest);
    %load('AR-50.mat')
    [Arows, Acols] = size(A);
    %degfree = Arows - Acols;
    %mse_pre = ((res)'*(res))/((.005/dt)^2*ARfac)/degfree
    %mse = (b-A*mest)'*pinv(AR)*A'*(b-A*mest)*ARfac
    mse = (b-A*mest)'*pinv(AR)*A'*(b-A*mest)  % simplifies because ARfac = 1 
    %clear AR
end
% irun = 1;
% 
% % find likelihood and prior for bayes factor
% if bi == 1
%     llikelihood_dT = -1/2*mse
%     lprior_dT = -1/2*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)
%     lpost_dT = lprior_dT*llikelihood_dT
% end
% 
% if bi == 2
%     llikelihood_dP = -1/2*mse
%     lprior_dP = -1/2*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)
%     lpost_dP = lprior_dP*llikelihood_dP
% end
% end

return



