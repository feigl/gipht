function varargout = funnull(varargin)
%function varargout = funnull(varargin)
%
%% return nulls

% number of arguments entered as input
nin = nargin;

% number of arguments to return as output
nout = nargout;

% decide what to do based on number of arguments

% --------------------- INITIALIZE 0 ---------
% call looks like this:
% [pnames,pscl] = feval(fitfun,me);
% [pnames,pscl] = feval(fitfun,me,datafilename);
% just return the names of the parameters in this model
if (nin == 1 || nin == 2) && nout == 2 && isstruct(varargin{1}) == 0
    fprintf(1,'Naming parameters in fitting function fitfun = %s\n',mfilename);
    
    pnames = {'Null'};
    pscl   = NaN;
    
    varargout(1) = {pnames};
    varargout(2) = {pscl};
    
    return;
    % --------------------- INITIALIZE 1 -------------
    % return with empty range vector
    %  call looks like this:
    % [rng,TST] = feval(fitfun,DST,PST);
elseif nin == 2 && nout == 2 && isstruct(varargin{1}) == 1 && isstruct(varargin{2}) == 1
    % and temporary storage structure TST
    fprintf(1,'Resetting TST in fitting function fitfun = %s\n',mfilename);
    
    % data structure
    DST = varargin{1};
    % parameter structure
    PST = varargin{2};
    
    
    % assume that we have only one type of data
    idatatype1 = DST.idatatype(1);
    
    % construct temporary storage structure
    TST.idatatype1 = idatatype1;
    
    
    % pass storage back to calling routine for future use
    varargout(1) = {rng};
    varargout(2) = {TST};
    return
    % --------------------- EVALUATE THE FITTING FUNCTION -------------
    % call looks like this:
    %   rng            = feval(fitfun,DST,PST,TST); % return set of scalar ranges
    %   [rng,DST]      = feval(fitfun,DST,PST,TST); % also return DST structure
    % containing modeled vector [DST.mx, DST.my, DST.mz]
    %elseif nin == 3 && nout == 1
elseif nin == 3 ...
        && (nout == 1 || nout == 2) ...
        && isstruct(varargin{1}) == 1 ...
        && isstruct(varargin{2}) == 1 ...
        && isstruct(varargin{3}) == 1
    
    % disp DD1; size(DD)
    % data structure
    DST = varargin{1};
    % parameter structure
    PST = varargin{2};
    % temporary storage structure
    TST = varargin{3};
    
    ndata = numel(DST.i);
    
    
    
    % return the modeled values
    rng0 = zeros(ndata,1);
    
    switch nout
        case 1
            varargout(1) = {rng0};
            % displacement vector if requested
        case 2
            varargout(1) = {rng0};
            DST.phamod = colvec(rng0);
            DST.mx = zeros(ndata,1);
            DST.my = zeros(ndata,1);
            DST.mz = zeros(ndata,1);
            varargout(2) = {DST};
        otherwise
            error(sprintf('Unknown value of nout = %d\n',nout));
    end
else
    error(sprintf('ERROR: nin = %d and nout = %d fitting function named %s\n'...
        ,nin,nout,mfilename));
end

return



