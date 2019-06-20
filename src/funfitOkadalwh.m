function varargout = funfitOkadalwh(varargin)
%function varargout = funfitOkada3(varargin)
%
%% separable model for range change
%   p       == parameter vector (1 column)
%   xyzm   == easting, northing, up in m [3 rows]
%   tepochs == unique times in decimal years [1 row]
%   bpest  == orbital separation w.r.t. virutal orbit in meters [1 row]
%   dops    == doppler separation w.r.t. virtual orbit in PRF [1 row]
%   DD      == differencing matrix of 1s and 0s
%   unitv   == unit vector pointing from ground to sat: E, N, U (dimless) [3 elements]
%   ippix1  == array of pixel pointers [np x 1]
%
% return range field rng which has ndata fields
% This version uses pd to make a per-epoch contribution
% This version includx 2 Mogi sources and 2 Okada sources, and differences
% the two Okada fields to simulate poroelastic effects
%
% Kurt Feigl 2008-JUL-24 optimize for speed
% 2009-MAR-28 handle arbitrary locations
% 2009-APR-02 calculate gradients with respect to mean coordinate
% 2009-JUN-18 rng in DN such 256 DN = 1 cycle
% 2009-JUL-15 keep partial
% 2009-OCT-21 skip ipass
% 2010-JAN-11 use east gradient in radians per pixel if ifast = 100 or 101
% 2010-OCT-05 use structures on input and output
% 2010-NOV-11 add PST.scale
% 2011-JUN    add Yang Model in place of Okada2 (Aaron)
% 2011-JUL    add penny crack
% 2011-SEP    add Poroelastic by subtracting 2 Okadas
% 2012-JAN-11 use CENTRAL difference for gradients in Mogi model
% 2012-JAN-14 use CENTRAL difference for gradients in all models
% 2012-JUN-25 return displacement vector, too
% 2012-NOV-26 call comsol
% 2014-MAY-14 abandon Abaqus, use general 3-D Comsol
% 2015-APR-19 replace Cervelli disloc routine with Beauducel okada85
% 2017-JUL-24 fix bug found by Helene - Kurt
% 2017-AUG-02 adapt for volumetric strain in 3 orthogonal dikes

% fprintf(1,'Entering fitting function fitfun = %s\n',mfilename);

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
    
    % number of epochs
    me = varargin{1};
    
    if nin == 2
        datafilename = varargin{2};
    else
        datafilename = '';
    end
    
    % start list of parameter names
    [pnames,pscl]=get_parameter_names_for_epochwise(me);
    mp = numel(pnames);  
    % continue list of parameter names
    j = mp;
    j=j+1;pnames{j} = sprintf('Reference Epoch in years       '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Poisson Ratio dimless          '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Volume Change in cubic meters  '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Easting in m            '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Northing in m           '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Elevation in m          '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 dX in m                 '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 dY in m                 '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 dZ in m                 '); pscl(j)=1.0E-3;
    

    % for debugging
    %pscl = ones(size(pscl));
    
    pnames = truncate_parameter_names(pnames);
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
    
    % number of data
    ndata = numel(DST.phaobs);
    
    % initialize
    rng1 = zeros(ndata,1);
    rng  = zeros(ndata,1);
    
    % number of parameters
    % mparam = numel(PST.p1); % number of parameters
    mparam = PST.mparam;
    if mparam ~= numel(PST.p1)
        error('Parameter miscount.')
    end
    
    % number of epochs
    % 20130624 me = numel(get_parameter_index('epoch',PST.names))/11;
    % WARNING the number of parameters per epoch is SET HERE 
    parameters_per_epoch = 11;
    me = numel(get_parameter_index('epoch',PST.names))/parameters_per_epoch;
    if round(me) ~= me
        error(sprintf('Number of epochs is NOT INTEGER %f\n',me));
    end
    
    % number of pairs
    np = numel(unique(DST.k));
    
    % index of first pair
    kpair = DST.k(1);
    
    % select row corresponding to this pair
    %DD1 = DD(kpair,:);
    
    DD  = zeros(ndata,me);
    DDM = zeros(ndata,me); % for master
    DDS = zeros(ndata,me); % for slave
    
%     fprintf(1,'Starting to build DD matrices for master and slave\nLine: ');
    for i = 1:ndata
       %fprintf(1,'%d ',i);
        % Row corresponding to this pair
        DD( i,DST.kmast(i)) = -1; % subtract master
       %DDS(i,DST.kmast(i)) = -1; % subtract master 20130629 - BUG
        DDM(i,DST.kmast(i)) = -1; % subtract master 20130629 - REPAIR
        DD( i,DST.kslav(i)) = +1; % from slave
        DDS(i,DST.kslav(i)) = +1; % from slave
        %DD(i,:)
    end
%     fprintf(1,'\nDone building DD matrices for master and slave\n');

%     disp 'First row of DD matrix';
%     DD(1,:)
%     disp 'First row of DDM matrix';
%     DDM(1,:)
%     disp 'First row of DDS matrix';
%     DDS(1,:)
    
    %     % remember this pair
    %     kmast = unique(DST.kmast);
    %     kslav = unique(DST.kslav);
    
    % mean location coordinates in meters
    %     xmean = mean(DST.x);
    %     ymean = mean(DST.y);
    %     zmean = mean(DST.z);
    
    % assume that we have only one type of data
    %fprintf(1,'Starting to check data types\n');

    idatatype1 = DST.idatatype(1);
    for i=1:ndata
        if DST.idatatype(i) ~= idatatype1
            fprintf(1,'idatatype differs for i = %d %d %d\n'...
                ,i,idatatype1,DST.idatatype(i));
        end
    end
%     fprintf(1,'Done checking data types\n');
    %rng = zeros(1,ndata);
    
    % constant for scaling baseline
    %bscale = -2 / 830.e3/ sind(23); % for ERS to yield m
    % 2011-JUN-24
    %incidence angle in radians from vertical
    bscale = acos(DST.uvz ./((DST.uvx).^2+(DST.uvy).^2+(DST.uvz).^2));
    %     fprintf(1,'Bscale (incidence angle in degrees) min, max, mean: %10.4f %10.4f %10.4f\n'...
    %         ,180.*min(bscale)/pi,180.*max(bscale)/pi,180.*mean(bscale)/pi);
    %disp 'bscale'; size(bscale)
    bscale = bscale - mean(bscale);
    
     
%     fprintf(1,'Starting to build TST structure\n');

    % construct temporary storage structure
    TST.idatatype1 = idatatype1;
    %TST.idatatype  = idatatype;
    TST.me         = me;
    TST.np         = np;
    TST.mparam     = mparam;
    TST.ndata      = ndata;
    TST.bscale     = bscale;
    TST.dd         = DD;
    TST.ddm        = DDM;
    TST.dds        = DDS;
    TST.p0         = PST.p0;
    TST.p1         = PST.p1;

%     fprintf(1,'Done building TST structure\n');

    
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

    
        %fprintf(1,'Evaluating fitting function fitfun = %s\n',mfilename);

%         fprintf(1,'Reading varargin\n');

    % disp DD1; size(DD)
    % data structure
    DST = varargin{1};
    % parameter structure
    PST = varargin{2};
    % temporary storage structure
    TST = varargin{3};

%     fprintf(1,'Unpacking TST\n');

    %% convert string to function handle
    timefun = str2func(PST.timefun);


    % unpack storage structure - must match above
    idatatype1 = TST.idatatype1;
    %idatatype  = TST.idatatype;
    me         = TST.me;
    np         = TST.np;
    mparam     = TST.mparam;
    ndata      = TST.ndata;
    bscale     = TST.bscale;
    DD         = TST.dd;
    DDM        = TST.ddm;
    DDS        = TST.dds;
    
%     fprintf(1,'Unpacking parameter vector\n');

    % get values of parameters
    p = colvec(PST.p1);
    pt =colvec(p( 0*me+1: 0*me+me)); % epoch
    px =colvec(p( 1*me+1: 1*me+me));
    py =colvec(p( 2*me+1: 2*me+me));
    pz =colvec(p( 3*me+1: 3*me+me));
    poh=colvec(p( 4*me+1: 4*me+me)); % parameters describing orbit position - horizontal  component
    poa=colvec(p( 5*me+1: 5*me+me)); % parameters describing orbit position - along-track component
    pov=colvec(p( 6*me+1: 6*me+me)); % parameters describing orbit position - vertical    component
    pvh=colvec(p( 7*me+1: 7*me+me)); % parameters describing orbit velocity - horizontal  component
    pva=colvec(p( 8*me+1: 8*me+me)); % parameters describing orbit velocity - along-track component
    pvv=colvec(p( 9*me+1: 9*me+me)); % parameters describing orbit velocity - vertical    component
    pd =colvec(p(10*me+1:10*me+me)); % additive offset
    pg =colvec(p(11*me+1:  mparam)); % geophysical parameters
 
    %     for ikp = 1:numel(pg)
    %         fprintf(1,'%d %s %12.4e\n',ikp,PST.names{ikp+koffset},pg(ikp));
    %     end
    
    %% pointer index into array pg
    kp = 0;

    %     j=j+1;pnames{j} = sprintf('Reference Epoch in years       '); pscl(j)=1.0;
    kp=kp+1;tquake = pg(kp);
    %     j=j+1;pnames{j} = sprintf('Poisson Ratio dimless          '); pscl(j)=1.0E-3;
    kp=kp+1;nu = pg(kp);
    %   j=j+1;pnames{j} = sprintf('Volume Change in cubic meters  '); pscl(j)=1.0E-3;
    kp=kp+1;F.dV = pg(kp);   
    %     j=j+1;pnames{j} = sprintf('Okada3 Easting in m            '); pscl(j)=1.0E-3;
    kp=kp+1;F.x = pg(kp);
    %     j=j+1;pnames{j} = sprintf('Okada3 Northing in m           '); pscl(j)=1.0E-3;
    kp=kp+1;F.y = pg(kp);
    %     j=j+1;pnames{j} = sprintf('Okada3 Elevation in m          '); pscl(j)=1.0E-3;
    kp=kp+1;F.z = pg(kp);
    %     j=j+1;pnames{j} = sprintf('Okada3 dX in m                 '); pscl(j)=1.0E-3;
    % all 3 dimensions are different
    kp=kp+1;F.dx = pg(kp);
    kp=kp+1;F.dy = pg(kp);
    kp=kp+1;F.dz = pg(kp);
     
    
    % check dimensions
    [nrDD,ncDD] = size(DD);
    if  nrDD ~= ndata || ncDD ~= me ...
            || numel(px) ~= me || numel(py) ~= me || numel(pz) ~= me || numel(pd) ~= me  ...
         
        me
        ndata
        disp DD; size(DD)
        disp px; size(px)
        disp py; size(py)
        disp pz; size(pz)
        disp pt; size(pt)
        disp pd; size(pd)
        error('Dimension error.');
    end
    
    %% nuisance parameters first
    if idatatype1 == -1
        % east component of gradient only
        nui = (DD * px) .* DST.dx;
    else
        %        nuisance parameters 3 components of gradient plus offset
        nui = (DD * px) .* DST.x0 ...
            + (DD * py) .* DST.y0 ...
            + (DD * pz) .* DST.z0 ...
            + (DD * pd) .* DST.mpercy; % additive constant
    end
    
    %% baseline term for orbits
    bas = (DDM * poh) .* DST.orbm1 ...
        + (DDM * poa) .* DST.orbm2 ...
        + (DDM * pov) .* DST.orbm3 ...
        + (DDS * poh) .* DST.orbs1 ...
        + (DDS * poa) .* DST.orbs2 ...
        + (DDS * pov) .* DST.orbs3 ...
        + (DDM * pvh) .* DST.orbm4 ...
        + (DDM * pva) .* DST.orbm5 ...
        + (DDM * pvv) .* DST.orbm6 ...
        + (DDS * pvh) .* DST.orbs4 ...
        + (DDS * pva) .* DST.orbs5 ...
        + (DDS * pvv) .* DST.orbs6;
    
    
   %% call the Okada routine 
    switch idatatype1
        case {0,2}  % observable is phase or range
%            uokada = okada85_wrapper3(rowvec(vecx),rowvec(vecy),F(i),nu); % use Kurt's function
           uokada = okada85_wrapperlwh(rowvec(DST.x),rowvec(DST.y),F,nu); % use Kurt's function
        %dr = -1*(unitv_east * uENZ(1,:) + unitv_north * uENZ(2,:) + unitv_up * uENZ(3,:));
            
        %case -1  % observable is gradient
        otherwise
            warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
            uokada = zeros(3,ndata);
    end
    

    %% ************************************
    % DONE WITH INDIVIDUAL PARTS OF DISPLACEMENT FIELD
    % ************************************
    % sum all RATES for separable part
    u = uokada;
%     disp 'u'; size u
    
    % position-dependent field of range change rate values in m/yr
    gmod = -1.0* (colvec(u(1,:)).*DST.uvx ...
                + colvec(u(2,:)).*DST.uvy ...
                + colvec(u(3,:)).*DST.uvz);
 
    %time function in a column vector
    %tdif = DD * colvec(pt);
    %tdif = DD * time_function(pt, tquake);
    %fprintf(1,'Time difference in years %10.4f\n',tdif);
    tdif = DD * timefun(pt, tquake);
    %fprintf(1,'Time difference in years %10.4f\n',tdif);

    
    %% combine time and space dependence in a column vector
    rng0 = tdif .* gmod;  % modeled range change in m
    
    %% add nuisance contribution
    rng1 = rng0 + nui + bas ;  

    %% create observable quantity 
    switch idatatype1
        case 0  % observable is phase in radians
            rng0 = 2.0 * pi * rng0 ./ DST.mpercy;
            rng =  2.0 * pi * rng1 ./ DST.mpercy;
%        case -1  % observable is dimensionless gradient
%             rng =  rng1;
        case 2  % observable is range change in meters
            rng = rng1;
        otherwise
            warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
            rng = zeros(ndata,1);
    end
 
    %% print statements for debugging
%     % replace NaN with zeros?
%     inan = find(isfinite(rng0) == 0);
%     rng0(inan) = 0;
    
%    fprintf(1,'min = %12.4f max = %12.4f \n',min(min(rng1)),max(max(rng1)));
    
    
%     fprintf(1,'Minimal values in meters for xyzm %10.1f %10.1f %10.1f \n',min(DST.x),min(DST.y),min(DST.z));
%     fprintf(1,'Maximal values in meters for xyzm %10.1f %10.1f %10.1f \n',max(DST.x),max(DST.y),max(DST.z));
% % %     fprintf(1,'Minimal value  in meters for rRad_Vel %10.4f\n',min(rRad_Vel));
% % %     fprintf(1,'Maximal value  in meters for rRad_Vel %10.4f\n',max(rRad_Vel));
%     fprintf(1,'Minimal value  in meters for rng  %10.4f\n',min(rng));
%     fprintf(1,'Maximal value  in meters for rng  %10.4f\n',max(rng));
%     fprintf(1,'Minimal value  in meters for bas  %10.4f\n',min(bas));
%     fprintf(1,'Maximal value  in meters for bas  %10.4f\n',max(bas));
%     fprintf(1,'i    Time Function pt(i)\n');
%     for i = 1:numel(pt)
%         fprintf(1,'%5d %12.6f\n',i,pt(i));
%     end
    

%     fprintf(1,'Done calculating. Returning...\n');

    %% return the modeled values
    switch nout
        case 1
            varargout(1) = {rng};
            % displacement vector if requested
        case 2
            varargout(1) = {rng0};
            DST.phamod = colvec(rng);
            % neglect nuisance effects
            DST.mx = tdif.*colvec(u(1,:));
            DST.my = tdif.*colvec(u(2,:));
            DST.mz = tdif.*colvec(u(3,:));
            %         % include nuisance effects
            %         DST.mx = tdif.*colvec(u(1,:)) + (nui + bas) .* DST.uvx;
            %         DST.my = tdif.*colvec(u(2,:)) + (nui + bas) .* DST.uvy;
            %         DST.mz = tdif.*colvec(u(3,:)) + (nui + bas) .* DST.uvz;
            varargout(2) = {DST};
        otherwise
            error(sprintf('Unknown value of nout = %d\n',nout));
    end
else
    error(sprintf('ERROR: nin = %d and nout = %d fitting function named %s\n'...
        ,nin,nout,mfilename));
end

return


