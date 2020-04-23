function [ x, c, s] = pinveb(A, tol, index)
% [ x, c, s] = pinveb(A, tol, index)
% Pseudoinverse with condion number 
% Finds the pseudoinverse of A using SVD, with option for truncation  
% INPUTS:  
%   A   - matrix
%   tol - tolerance level for singular values or cut-off value for
%         number of singular values; in the latter case, enter 'index'
%         after tol
%   p   - cut-off value for number of singular values 
%   Leave tol and p unspecified for automatic truncation based on built in tol
% 
% OUTPUTS:  
%   x - pseudoinverse of A 
%   c - condition number of pseudoinverse 
%
% This function also returns a plot of the singular value spectrum with 
% truncation index identified
%
% Elena C. Baluyut, UW-Madison
%2015 - 03 - 08

% Compute SVD
[U S V] = svd(A);

% Determine number of singluar values
[mtest ntest] = size(S);
if ntest == 1
    s = S;
%    disp('test 1')
else
    s = diag(S); % store singular values in a vector
end

n = 1:numel(s); % store indices of singluar values

% Assign tolerance and cut-off values 
if nargin < 2 
    % if tolerance isn't specified, use default
    tol = max(size(A)) * eps(norm(s,inf)); 
    sp = s(s > tol);
    p = numel(sp);
elseif nargin == 2 
    % find cut-off based on tolerance level
    sp = s(s > tol);
    p = numel(sp);
elseif nargin == 3 
    %cut-off specified, select subset
    sp = s(1:tol);
    p = numel(sp);
else
    error('too many inputs')
end

pn = 1:p; % set indices for truncated singular values

% Truncate SVD 
Up = U(:, pn);
Sp = S(pn, pn);
Vp = V(:, pn);

% Calculate pseudoinverse and condition number
x = Vp*inv(Sp)*Up';
c = cond(x);

% Plot singular value spectrum with truncation cut-off
% 20200205 Kurt Feigl - only plot if more than one singular value
if n > 1 
    figure
    plot(n, s/s(1), 'ko')
    hold on
    plot(p*ones(1, 10), linspace(min(s./s(1)), max(s./s(1)), 10), 'b--')
    hold off
    xlabel('i')
    ylabel('(s/s(1))')
    title('L curve for singluar values showing current truncation')
    axis tight
end

end

