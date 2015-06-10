function rng = funtaylor4(DST,PST,TST)
%function rng = funtaylor4(DST,PST,TST)
% approximate range using fourth order Taylor expansion
%  phamod=phamod0+matmul(G1,p-x0)+(1.0/2.0)*matmul(G2,(p-x0)**2)+(1.0/6.0)*matmul(G3,(p-x0)**3)+(1.0/24.0)*matmul(G4,(p-x0)**4)


% number of arguments entered as input
nin = nargin;

% number of arguments to return as output
nout = nargout;

if nargin == 3 && nargout == 1
    % get dimensions
    mparam     = PST.mparam;
    ndata      = numel(DST.i);
    [nrows, mcols] = size(TST.partial_wrt_1param);
    
    if nrows ~= ndata
        error(sprintf('Number of rows (%d) DOES NOT EQUAL number of data points (%d)\n',nrows,ndata));
    end
    if mcols ~= mparam
        error(sprintf('Number of columns (%d) DOES NOT EQUAL number of parameters (%d)\n',mcols,mparam));
    end
    
%     nnan = numel(find(isfinite(TST.partial_wrt_1param) == 0));
%     if nnan > 0
%         error('Found NaN values in matrix of partial derivatives');
%     end
%     
    % calculate (unwrapped) range change in radians
%   rng = DST.phamod + TST.partial_wrt_1param * ((PST.p1 - PST.p0) ./ PST.scale);
    rng = DST.phamod + TST.partial_wrt_1param * (PST.p1 - PST.p0);
    
    % return a ROW vector
    %rng = rowvec(rng);
    %size(rng)
else
    error(sprintf('ERROR in nin = %d and nout = %d fitting function named %s\n'...
        ,nin,nout,mfilename));
end
       
return

