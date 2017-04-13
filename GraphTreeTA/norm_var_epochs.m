function [ Wd, Wp ] = norm_var_epochs(Q, trees, dsig)
% function [ Wp ] = norm_inc_matrix(Q, trees, dsig )
% Computes relative variances of epoch-wise measurements by regularizing
% the edge-vertex incidence matrix with a constraint of zero-mean per
% component.
%
% INPUTS:
%   Q - edge-vertex incidence matrix
%   trees - matrix of epoch indices per component in the data set
%   dsig - vector of uncertainties in pair-wise measurements
%   nrm - 1 for normalizing, 0 for epoch-wise misfit estimation for MST
%
% OUTPUT:
%   Wp - vector of relative variances of epoch-wise measurements
%
% Elena C. Baluyut, UW-Madison
% 2015-07-26

% Initialize

Qp = Q;
[ndat, mdummy] = size(Q);
dsigp = dsig;

[ntrees, nepochs] = size(trees);
% if nrm == 0
%     Vp = incidence_to_cov(Qp, dsigp);
%     Wp = pinv(Qp'*Qp)*Qp'*Vp*Qp*pinv(Qp'*Qp);
% else
    Qp(end+1:end+ntrees, :) = 0;
    dsigp(end+1:end+ntrees) = 1;
    
    % Loop through trees
    
    for n = 1:ntrees
        iok = isfinite(trees(n,:)); % locate epochs in component
        nodes = trees(n,iok);   % identify indices of epochs in component
        Qp(ndat+n, nodes) = 1/numel(nodes); % assign extra row with averaging constraint
        % check for accuracy
        if isequal(find(Qp(ndat+n, :)), nodes) == 0
            error('incorrect assignment of averaging elements');
        end
    end
    
    % Define new data Covariance matrix
    
    Vp = incidence_to_cov(Qp, dsigp);
    Wp = inv(Qp'*Qp)*Qp'*Vp*Qp*inv(Qp'*Qp); % covariance of epoch-wise measurements
% end
Wd = diag(Wp);

end

