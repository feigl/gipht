function [ pest, psig, mse, Vx, index ] = lsqr_tikh( G, V, d, L, beta_range )
% function [ pest, psig, mse, Vx, index, minnorm ] = lsqr_tikh( G, V, tu, d, L, beta_range )
% Performs interative weighted least squares inversion for Tikhonov
% regularization
%
% INPUTS:
%   G -            design matrix (size n-by-m)
%   V -            covariance matrix of pair-wise data (size n-by-n)
%   d -            vector of pair-wise data (size n-by-1)
%   L -            Tikhonov regularization matrix (size (m, m-1, or m-2)-by-m,
%                  depending on order (m for order 0, m-1 for order 1, m-2 for order 2)
%   beta_range -  range for interation for beta (regularization parameter)
%
% OUTPUTS:
%   pest -        least-squares estimated solution for unknown parameter vector, size
%                 [m - by 1]
%   psig -        vector of uncertainties for estimated solution
%   mse -         mean squared error (or variance of unit weight) of estimated
%                 solution
%   Vx -          covariance matrix of estimated parameters
%   index -       index of regularization parameter used from alph
%    
% Elena C. Baluyut, UW-Madison
% 2015-07-19

% Initialize variables
[ndat, mparams] = size(G);

% Define storage matrices
M = zeros(mparams, numel(beta_range)); % used to store estimated parameters
rho = zeros(size(beta_range));         % used to store norm of residual
eta = zeros(size(beta_range));         % used to store norm of Tikhonov fit

% Iterate through alphas
for i = 1:numel(beta_range)
bL = beta_range(i)^2*L;   % define beta*L
[nl ml] = size(bL);       % find size of beta*L
Ze = zeros(nl, 1);        % set corresponding zero vector

% define new matrices for system
GL = [sqrtm(pinv(V))*G; bL];
dL = [sqrtm(pinv(V))*d; Ze];

%perform inversion
m = pinv(GL'*GL)*GL'*dL;
M(:, i) = m;                % store estimated parameters
rho(i) = norm(G*m-d);       % store norm of residual
eta(i) = norm(L*m);         % store norm of Tikhonov fit

end

% Iterate through plots 
sats =0;
while sats < 1
    
    % Plot L curve
    clf
    plot(rho, eta, 'ko')
    xlabel('norm(G*m-d)')
    ylabel('norm(L*m)')
    hold on
    
    % FOR INTERACTIVE DECISION, REMOVE COMMENT ON NEXT LINE
    %interactive = 1;
    interactive = 0;
    if interactive == 1
        % User defines index for parameterization
        %index = input('choose index of alpha: ')
        
        % User desides to keep index or choose a new one
        check = input('continue with this alpha? [y/n]: ', 's')
        if strcmp(check, 'y') == 0 && strcmp(check, 'n') == 0
            disp('incorrect input, make sure not to add spaces')
        end
       
    else
        % FOR OKMOK EXAMPLE ONLY
        index = 2
        check = 'y'     
    end
    
    plot(rho(index), eta(index), 'r*') % plot chosen index on L curve
    hold off 
    

     sats = strcmp(check, 'y'); % repeats process until index is satisfactory, 
end

% Compute statistics similarly to other parameterizations
pest = M(:, index);
if isreal(pest) ~= 1
    warning('pest has an imaginary part with magnitude %f Taking real part',nanmax(imag(pest)));
    pest = real(pest);
end

mse = d'*(pinv(V) - pinv(V)*G*pinv(G'*pinv(V)*G)*G'*pinv(V))*d./(ndat-mparams);

Vx = pinv(G'*pinv(V)*G)*mse;

% standard deviation of parameters
psig = sqrt(diag(pinv(G'*pinv(V)*G)*mse));
if isreal(psig) ~= 1
    warning('psig has an imaginary part with magnitude %f Taking real part',nanmax(imag(psig)));
    psig = real(psig);
end

return



