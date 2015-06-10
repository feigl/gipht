function varargout = funseparable27(varargin)
%function varargout = funseparable27(varargin)
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

% number of arguments entered as input
nin = nargin;

% number of arguments to return as output
nout = nargout;

% decide what to do based on number of arguments

% --------------------- INITIALIZE 0 ---------
% call looks like this:
% [pnames,pscl] = feval(fitfun,me);
% just return the names of the parameters in this model
if nin == 1 && nout == 2 && isstruct(varargin{1}) == 0
    fprintf(1,'Naming parameters in fitting function fitfun = %s\n',mfilename);
    
    % number of epochs
    me = varargin{1};
    
    % start list of parameter names
    [pnames,pscl]=get_parameter_names_for_epochwise(me);
    mp = numel(pnames);
    
    % continue list of parameter names
    j = mp;
    j=j+1;pnames{j} = sprintf('Mogi1 Easting in m             '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Mogi1 Northing in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Mogi1 Depth in m               '); pscl(j)=1.0E4;
    j=j+1;pnames{j} = sprintf('Mogi1 Volume Increase in m3    '); pscl(j)=1.0E6;
    j=j+1;pnames{j} = sprintf('Mogi2 Easting in m             '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Mogi2 Northing in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Mogi2 Depth in m               '); pscl(j)=1.0E4;
    j=j+1;pnames{j} = sprintf('Mogi2 Volume Increase in m3    '); pscl(j)=1.0E6;
    j=j+1;pnames{j} = sprintf('Okada1 Length in m             '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada1 Width in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada1 Depth in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada1 Negative Dip in deg     '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada1 Strike CCW from N in deg'); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada1 Easting in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada1 Northing in m           '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada1 RL Strike Slip in m     '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada1 Downdip Slip in m       '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada1 Tensile Opening in m    '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada2 Length in m             '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada2 Width in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada2 Depth in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada2 Negative Dip in deg     '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada2 Strike CCW from N in deg'); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada2 Easting in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada2 Northing in m           '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada2 RL Strike Slip in m     '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada2 Downdip Slip in m       '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada2 Tensile Opening in m    '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('YangPS Easting in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('YangPS Northing in m           '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('YangPS Depth in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('YangPS Excess Pressure in Pa   '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('YangPS semimajor axis a  m     '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('YangPS semiminor axis b in m   '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('YangPS Azimuth deg CCW from N  '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('YangPS Plunge in degrees       '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('YangPS Empty1                  '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('YangPS Empty2                  '); pscl(j)=1.0E-3;
    % the following parameters are unchanged
    j=j+1;pnames{j} = sprintf('Poisson Ratio dimless          '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Shear Modulus in Pa            '); pscl(j)=1.0E10;
    j=j+1;pnames{j} = sprintf('Reference Epoch in years       '); pscl(j)=1.0;
    %j=j+1;pnames{j} = sprintf('Young s Modulus in Pa          '); pscl(j)=1.0E9;
    j=j+1;pnames{j} = sprintf('Young s Modulus in Pa for Pinel'); pscl(j)=1.0E9;
    %j=j+1;pnames{j} = sprintf('Density of ice in kg per m3    '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Density of load in kg per m3   '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Grav Acc in m per s per s      '); pscl(j)=1.0;
    % new 2011-JUL-23
    j=j+1;pnames{j} = sprintf('Poisson Ratio dimless for Pinel'); pscl(j)=1.0E10;
    j=j+1;pnames{j} = sprintf('SunDisk1 Easting in m          '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk1 Northing in m         '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk1 Depth in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk1 Excess Pressure in Pa '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('SunDisk1 Radius a in m         '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk2 Easting in m          '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk2 Northing in m         '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk2 Depth in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('SunDisk2 Excess Pressure in Pa '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('SunDisk2 Radius a in m         '); pscl(j)=1.0E3;
    % Added 2011-SEP-21
    j=j+1;pnames{j} = sprintf('Poisson Ratio dimless drained  '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Length in m             '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada3 Width in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada3 Depth in m              '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada3 Negative Dip in deg     '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada3 Strike CCW from N in deg'); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('Okada3 Easting in m            '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada3 Northing in m           '); pscl(j)=1.0E3;
    j=j+1;pnames{j} = sprintf('Okada3 RL Strike Slip in m     '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Downdip Slip in m       '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Okada3 Tensile Opening in m    '); pscl(j)=1.0E-3;
    % Added 2012-NOV-26
    j=j+1;pnames{j} = sprintf('CS Easting in m                '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('CS Northing in m               '); pscl(j)=1.0e3;
%     j=j+1;pnames{j} = sprintf('CS POR Porosity                '); pscl(j)=1.0;
%     j=j+1;pnames{j} = sprintf('CS PER Permeability in m2      '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('CS POR Log10 of Porosity       '); pscl(j)=1.0;
    j=j+1;pnames{j} = sprintf('CS PER Log10 Permeability m2   '); pscl(j)=1.0;
    % Added 2013-05-21
    j=j+1;pnames{j} = sprintf('AQ Excess Pressure in Pa       '); pscl(j)=1.0e6;
    j=j+1;pnames{j} = sprintf('AQ Easting in m                '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('AQ Northing in m               '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('AQ Depth in m                  '); pscl(j)=1.0;
    % Added 2013-06-23
    j=j+1;pnames{j} = sprintf('Rad Vel Center Easting in m    '); pscl(j)=1.0e6;
    j=j+1;pnames{j} = sprintf('Rad Vel Center Northing in m   '); pscl(j)=1.0e6;
    j=j+1;pnames{j} = sprintf('Rad Vel Evel in m per yr       '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('Rad Vel Nvel in m per yr       '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('Rad Vel Uvel in m per yr       '); pscl(j)=1.0e3;
    j=j+1;pnames{j} = sprintf('Rad Vel dVdR in inverse yr     '); pscl(j)=1.0e3;
    
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
    %fprintf(1,'Resetting TST in fitting function fitfun = %s\n',mfilename);
    
    % data structure
    DST = varargin{1};
    % parameter structure
    PST = varargin{2};
    
    % number of epochs
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
    for i = 1:ndata
        % Row corresponding to this pair
        DD( i,DST.kmast(i)) = -1; % subtract master
       %DDS(i,DST.kmast(i)) = -1; % subtract master 20130629 - BUG
        DDM(i,DST.kmast(i)) = -1; % subtract master 20130629 - REPAIR
        DD( i,DST.kslav(i)) = +1; % from slave
        DDS(i,DST.kslav(i)) = +1; % from slave
        %DD(i,:)
    end
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
    idatatype1 = DST.idatatype(1);
    for i=1:ndata
        if DST.idatatype(i) ~= idatatype1
            fprintf(1,'idatatype differs for i = %d %d %d\n'...
                ,i,idatatype1,DST.idatatype(i));
        end
    end
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
    
    % do once to initialize
    upinel1a = zeros(3,ndata);
    upinel1b = zeros(3,ndata); 
    uabaqus1a = zeros(3,ndata);
    uabaqus1b = zeros(3,ndata);
    % Pointers for parameter names
    koffset   =  parameters_per_epoch * me;
    kAQ       = get_parameter_index('AQ_Excess_Pressure_in_Pa________',PST.names);
    kCS       = get_parameter_index('CS_Easting_in_m_________________',PST.names);
    kRad_Vel  = get_parameter_index('Rad_Vel_Center_Easting_in_m_____',PST.names);
    kYangPS   = get_parameter_index('YangPS_Easting_in_m_____________',PST.names);
    kSunDisk1 = get_parameter_index('SunDisk1_Easting_in_m___________',PST.names);
    kSunDisk2 = get_parameter_index('SunDisk2_Easting_in_m___________',PST.names);
    kPinel    = get_parameter_index('Young_s_Modulus_in_Pa_for_Pinel_',PST.names);
  
    % data file is named
    if numel(char(PST.datafilename)) > 0
        % data file exists
        if fexist(char(PST.datafilename)) == 1           
            % Need to initialize Pinel model
            if abs(PST.p0(get_parameter_index('Density_of_load_in_kg_per_m3____',PST.names))) > 0.
                [XI,YI,HI,mn] = pinel3d(char(PST.datafilename));
                switch idatatype1
                    case 0 % observable is phase
                        upinel1a = pinel3d(DST.x,DST.y,XI,YI,HI,mn);
                        upinel1b = zeros(3,ndata);
                    case -1 % observable is gradient: is the sign below flipped?
                        upinel1a = pinel3d(DST.x-DST.dx,DST.y,XI,YI,HI,mn);
                        upinel1b = pinel3d(DST.x+DST.dx,DST.y,XI,YI,HI,mn);
                    otherwise
                        warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                        upinel1a = zeros(3,ndata);
                        upinel1b = zeros(3,ndata);
                end
            end
            
            % Need to initialize Abaqus model
            %kAQ = get_parameter_index('AQ_Excess_Pressure_in_Pa________',PST.names);
            if abs(PST.p1(kAQ)) > 0.
                PSTAQ.mname    = PST.datafilename;
                PSTAQ.pressure = PST.p1(kAQ+0); % Pressure increase in Pascals
                PSTAQ.xcen     = PST.p1(kAQ+1); % Easting coordinate of Center in meters
                PSTAQ.ycen     = PST.p1(kAQ+2); % Northing coordinate of Center in meters
                PSTAQ.depth    = PST.p1(kAQ+3); % Depth of center in meters
                %PSTAQ.EYM = E;
                %PSTAQ.NUP = nu;             % Poisson's ratio
                %PSTAQ.verbose = 1;          % tell us what you are doing
                PSTAQ.verbose = 0;          % please be quiet
                              
                % Run the model from scratch             
                switch idatatype1
                    case 0 % observable is phase
                        uabaqus1a = abaqus2d_mogi1(PSTAQ,DST.x,DST.y);
                        uabaqus1b = zeros(3,ndata);
                    case -1 % observable is gradient
                        uabaqus1a = abaqus2d_mogi1(PSTAQ,DST.x-DST.dx,DST.y);
                        uabaqus1b = abaqus2d_mogi1(PSTAQ,DST.x+DST.dx,DST.y);
                    otherwise
                        warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                        uabaqus1a = zeros(3,ndata);
                        uabaqus1b = zeros(3,ndata);
                end
            end
        end
    end
    
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
    TST.upinel1a   = upinel1a;
    TST.upinel1b   = upinel1b;
    TST.uabaqus1a  = uabaqus1a;
    TST.uabaqus1b  = uabaqus1b;
    TST.parameters_per_epoch = parameters_per_epoch;
    TST.koffset    = koffset;
    TST.kAQ        = kAQ;
    TST.kCS        = kCS;
    TST.kYangPS    = kYangPS;
    TST.kSunDisk1  = kSunDisk1;
    TST.kSunDisk2  = kSunDisk2;
    TST.kRad_Vel   = kRad_Vel;
    TST.kPinel     = kPinel;
    
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
    upinel1a   = TST.upinel1a;
    upinel1b   = TST.upinel1b;
    uabaqus1a  = TST.uabaqus1a;
    uabaqus1b  = TST.uabaqus1b;
    koffset    = TST.koffset;
    kCS        = TST.kCS;
    kAQ        = TST.kAQ;
    kYangPS    = TST.kYangPS;
    kSunDisk1  = TST.kSunDisk1;
    kSunDisk2  = TST.kSunDisk2;
    kRad_Vel   = TST.kRad_Vel;
    kPinel     = TST.kPinel;    
    
    % get values of parameters
    p = colvec(PST.p1);
    pt =p( 0*me+1: 0*me+me); % epoch
    px =p( 1*me+1: 1*me+me);
    py =p( 2*me+1: 2*me+me);
    pz =p( 3*me+1: 3*me+me);
    poh=p( 4*me+1: 4*me+me); % parameters describing orbit position - horizontal  component
    poa=p( 5*me+1: 5*me+me); % parameters describing orbit position - along-track component
    pov=p( 6*me+1: 6*me+me); % parameters describing orbit position - vertical    component
    pvh=p( 7*me+1: 7*me+me); % parameters describing orbit velocity - horizontal  component
    pva=p( 8*me+1: 8*me+me); % parameters describing orbit velocity - along-track component
    pvv=p( 9*me+1: 9*me+me); % parameters describing orbit velocity - vertical    component
    pd =p(10*me+1:10*me+me); % additive offset
    pg =p(11*me+1:  mparam); % geophysical parameters
 
    %     for ikp = 1:numel(pg)
    %         fprintf(1,'%d %s %12.4e\n',ikp,PST.names{ikp+koffset},pg(ikp));
    %     end

    % pointer index into array pg
    %kp = 2*4 + 2*10 +  1;
    kp = 2*4 + 3*10 +  1;
    %  51 Poisson_Ratio_dimless___________
    nu = pg(kp);
    %  52 Shear_Modulus_in_Pa_____________
    kp=kp+1;
    mu = pg(kp);
    %  53 Reference_Epoch_in_years________
    kp=kp+1;
    tquake = pg(kp);
    %  54 Young_s_Modulus_in_Pa___________
    kp=kp+1;
    youngs = pg(kp);
    %  55 Density_of_load_in_kg_per_m3_____
    kp=kp+1;
    rho_ld = pg(kp);
    %  56 Grav._Acc._in_m_per_s_per_s_____
    kp=kp+1;
    gravity = pg(kp);
    %  57 Poisson_Ratio_for Pinel
    kp=kp+1;
    nu_pinel = pg(kp);
    %   Poisson_Ratio_dimless_drained___
    kp=56;
    nu2 = pg(kp);
    
    % check dimensions
    [ndum,mdum] = size(DD);
    if  ndum ~= ndata || mdum ~= me ...
            || numel(px) ~= me || numel(py) ~= me || numel(pz) ~= me || numel(pd) ~= me
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
    
    % nuisance parameters first
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
    
    % baseline term for orbits
    %bas = zeros(ndata,1);
%      poh
%      DDM(1,:) * poh + DDS(1,:) * poh
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
    
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % calculate Mogi displacments
    % Mogi    [u,e,t]=Mogi(volgeom, xloc, nu)
    %
    %  Computes surface displacements, strains, and tilts due to a Mogi source.
    %
    %  Inputs:
    %       volgeom = Mogi source geometry: East, North, Depth, Volume change
    %                 <length, length, length, volume>
    %          xloc = Matrix of local station coordinates <length>, stored columnwise,
    %                 i.e., east in first row, north in second row
    %            nu = Poisson's ratio
    %
    %  Output:
    %             u = matrix of displacements: Ux, Uy, Uz <length>
    %             e = matrix of strains: Exx, Exy, Eyy
    %             t = matrix of tilts: dUz/dx, dUz/dy
    %
    %  Notes: The term 'depth' denotes an unsigned length and should therefore always be
    %  given positive.  Keep your length units consistent! If you mix km and m you may get
    %  unexpected results, particularly with the strains and tilts.
    %
    %  05-17-98 Peter Cervelli.
    
    volgeom1(1) = pg(1);          % E coordinate of source in m
    volgeom1(2) = pg(2);          % N coordinate of source in m
    volgeom1(3) = pg(3);          % Depth of source in m
    volgeom1(4) = pg(4);          % Volume change of source in m^3
    volgeom2(1) = pg(5);          % E coordinate of source in m
    volgeom2(2) = pg(6);          % N coordinate of source in m
    volgeom2(3) = pg(7);          % Depth of source in m
    volgeom2(4) = pg(8);          % Volume change of source in m^3
    
    % call the function if volume is not zero
    switch idatatype1
        case 0  % observable is phase
            if abs(volgeom1(4)) > 0
                umogi1 = mogi(volgeom1, [DST.x';DST.y'], nu);
            else
                umogi1 = zeros(3,ndata);
            end
            
            if abs(volgeom2(4)) > 0
                umogi2 = mogi(volgeom2, [DST.x';DST.y'], nu);
            else
                umogi2 = zeros(3,ndata);
            end
        case -1 % observable is east gradient
            if abs(volgeom1(4)) > 0
                tl = mogi(volgeom1, [DST.x'-DST.dx';DST.y'], nu);
                tr = mogi(volgeom1, [DST.x'+DST.dx';DST.y'], nu);
                umogi1 = (tr - tl)/2.0;
            else
                umogi1 = zeros(3,ndata);
            end
            if abs(volgeom2(4)) > 0
                tl = mogi(volgeom2, [DST.x'-DST.dx';DST.y'], nu);
                tr = mogi(volgeom2, [DST.x'+DST.dx';DST.y'], nu);
                umogi2 = (tr - tl)/2.0;
            else
                umogi2 = zeros(3,ndata);
            end
        otherwise
            warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
            umogi1 = zeros(3,ndata);
            umogi2 = zeros(3,ndata);
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % CALCULATE OKADA DISLOCATIONS
    % OKADA MODEL PARAMETERS
    %
    % pg(5)   pg(6)    pg(7)  pg(8)  pg(9)  pg(10) pg(11)   pg(14)   pg(15)  pg(16)
    % pg(9)   pg(10)   pg(11) pg(12) pg(13) pg(14) pg(15)   pg(16)   pg(17)  pg(18)
    % pg(19)  pg(20)   pg(21) pg(22) pg(23) pg(24) pg(25)   pg(26)   pg(27)  pg(28)
    % Length - Width - Depth - Dip - Strike - East - North - Sslip - Dslip - Op
    % (km)     (km)    (km)    (deg) (deg)    (km)   (km)    (m)     (m)     (m)
    %                  upper   neg-
    %                  edge    ative
    
    % Call disloc only if slip is  non-zero
    % call the function if volume is not zero
    switch idatatype1
        case 0  % observable is phase
            kp = 9;
            % coseismic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                uokada1 = disloc(pg(kp:kp+9),[DST.x';DST.y'],nu);
            else
                uokada1 = zeros(3,ndata);
            end
            kp = 19;
            % coseismic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                uokada2 = disloc(pg(kp:kp+9),[DST.x';DST.y'],nu);
            else
                uokada2 = zeros(3,ndata);
            end
            kp = 57;
            % poroelastic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                uokada3 = disloc(pg(kp:kp+9),[DST.x';DST.y'],nu2 ) ...
                    -     disloc(pg(kp:kp+9),[DST.x';DST.y'],nu  );
            else
                uokada3 = zeros(3,ndata);
            end
            
        case -1  % observable is gradient
            kp = 9;
            % coseismic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                tl = disloc(pg(kp:kp+9),[DST.x'-DST.dx';DST.y'],nu);
                tr = disloc(pg(kp:kp+9),[DST.x'+DST.dx';DST.y'],nu);
                uokada1 = (tr - tl)/2.0;
            else
                uokada1 = zeros(3,ndata);
            end
            kp = 19;
            % coseismic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                tl = disloc(pg(kp:kp+9),[DST.x'-DST.dx';DST.y'],nu);
                tr = disloc(pg(kp:kp+9),[DST.x'+DST.dx';DST.y'],nu);
                uokada2 = (tr - tl)/2.0;
            else
                uokada2 = zeros(3,ndata);
            end
            kp = 57;
            % poroelastic
            if  sum(abs(pg(kp+7:kp+9))) > 1.e-3
                tl      = disloc(pg(kp:kp+9),[DST.x'-DST.dx';DST.y'],nu2 ) ...
                    -     disloc(pg(kp:kp+9),[DST.x'-DST.dx';DST.y'],nu  );
                tr      = disloc(pg(kp:kp+9),[DST.x'+DST.dx';DST.y'],nu2 ) ...
                    -     disloc(pg(kp:kp+9),[DST.x'+DST.dx';DST.y'],nu  );
                uokada3 = (tr - tl)/2.0;
            else
                uokada3 = zeros(3,ndata);
            end
        otherwise
            warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
            uokada1 = zeros(3,ndata);
            uokada2 = zeros(3,ndata);
            uokada3 = zeros(3,ndata);
    end
    
    % 20101215 AEMasters
    % Yang model for prolate spheroid
    %         xs    = as(1);     % center x
    %         ys    = as(2);     % center y
    %         z0    = as(3);     % center depth (positive)
    %         P     = as(4);     % excess pressure, mu*10^(-5) Pa
    %         a     = as(5);     % major axis, km
    %         b     = as(6);     % minor axis, km
    %         phi   = as(7);     % strike, rad  (0-2*pi)
    %         theta = as(8);     % plunge, rad  (0-pi)
    %     j=j+1;pnames{j} = sprintf('YangPS Easting in m            '); pscl(j)=1.0E3;
    %     j=j+1;pnames{j} = sprintf('YangPS Northing in m           '); pscl(j)=1.0E3;
    %     j=j+1;pnames{j} = sprintf('YangPS Depth in m              '); pscl(j)=1.0E3;
    %     j=j+1;pnames{j} = sprintf('YangPS Excess Pressure in Pa   '); pscl(j)=1.0E-3;
    %     j=j+1;pnames{j} = sprintf('YangPS semimajor axis a  m     '); pscl(j)=1.0E3;
    %     j=j+1;pnames{j} = sprintf('YangPS semiminor axis b in m   '); pscl(j)=1.0E3;
    %     j=j+1;pnames{j} = sprintf('YangPS Azimuth deg CCW from N  '); pscl(j)=1.0;
    %     j=j+1;pnames{j} = sprintf('YangPS Plunge in degrees       '); pscl(j)=1.0;
    %     j=j+1;pnames{j} = sprintf('YangPS Empty1                  '); pscl(j)=1.0E-3;
    %     j=j+1;pnames{j} = sprintf('YangPS Empty2                  '); pscl(j)=1.0E-3;
    
    %kp = 29;
    kp = kYangPS - koffset;
    if (abs(pg(kp+3)) > 1)
        matrl(1) = 2 * mu * nu / (1 - 2 * nu); %  1st Lame
        matrl(2) = mu;
        matrl(3) = nu;
        switch idatatype1
            case 0  % observable is phase
                [uYangPS(1,:), uYangPS(2,:), uYangPS(3,:)] = fcn_yangM(pg(kp:kp+7),DST.x',DST.y',matrl,DST.z');
            case -1 % observable is gradient
                [luY1, luY2, luY3] = fcn_yangM(pg(kp:kp+7),DST.x'-DST.dx',DST.y',matrl,DST.z');
                [ruY1, ruY2, ruY3] = fcn_yangM(pg(kp:kp+7),DST.x'+DST.dx',DST.y',matrl,DST.z');
                uYangPS(1,:) = (ruY1 - luY1)/2.0;
                uYangPS(2,:) = (ruY2 - luY2)/2.0;
                uYangPS(3,:) = (ruY3 - luY3)/2.0;
            otherwise
                warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                uYangPS = zeros(3,ndata);
        end
    else
        uYangPS = zeros(3,ndata);
    end
    
    % calculate displacements [east, north, up] at [xobs, yobs] due to a of height HI on grid XI, YI
    % Pinel et al. (2008) Geophys. J. Int. v. 169, pp 325-338
    % assume value of density (rho) is same everywhere
    % Kurt Feigl 2010-OCT-05
    %     gravity = 9.8;     % graviational acceleration in m/s/s
    %     rho_ld  = 0.9e3; % density of ice in kg/m^3
    %     youngs  = 60e9;    % Young's modulus in Pa
    %     loadfile = 'VaWEmism93-95.dat'
    %nu = 0.25;   % Poisson's ratio
    if rho_ld > 1.0
        switch idatatype1
            case 0  % observable is phase              
                upinel = pinel3d(gravity,youngs,nu_pinel,rho_ld,upinel1a);
            case -1 % observable is gradient
                tl=pinel3d(gravity,youngs,nu_pinel,rho_ld,upinel1a);
                tr=pinel3d(gravity,youngs,nu_pinel,rho_ld,upinel1b);
                upinel = (tr - tl)/2.0;
            otherwise
                warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                upinel = zeros(3,ndata);               
        end
    else
        upinel = zeros(3,ndata);
    end
    
    %SUN69   Deformation from penny-shaped crack in elastic half-space.
    %
    %  Author: Franois Beauducel <beauducel@ipgp.fr>
    %    Institut de Physique du Globe de Paris, 2009.
    %
    %  References:
    %   Sun, R. J. (1969). Theoretical size of hydraulically induced horizontal
    %     fractures and corresponding surface uplift in an idealized medium,
    %     J. Geophys. Res., 74, 5995-6011.
    %
    
    E = mu / (1 + 2.0*nu);  % Young's Modulus (un-numbered equation Sun p. 9997)
    
    %
    %
    %     0 j=j+1;pnames{j} = sprintf('SunDisk1 Easting in m          '); pscl(j)=1.0E3;
    %     1 j=j+1;pnames{j} = sprintf('SunDisk1 Northing in m         '); pscl(j)=1.0E3;
    %     2 j=j+1;pnames{j} = sprintf('SunDisk1 Depth in m            '); pscl(j)=1.0E3;
    %     3 j=j+1;pnames{j} = sprintf('SunDisk1 Excess Pressure in Pa '); pscl(j)=1.0E-3;
    %     4 j=j+1;pnames{j} = sprintf('SunDisk1 Radius a in m         '); pscl(j)=1.0E3;
    %     5 j=j+1;pnames{j} = sprintf('SunDisk2 Easting in m          '); pscl(j)=1.0E3;
    %     6 j=j+1;pnames{j} = sprintf('SunDisk2 Northing in m         '); pscl(j)=1.0E3;
    %     7 j=j+1;pnames{j} = sprintf('SunDisk2 Depth in m            '); pscl(j)=1.0E3;
    %     8 j=j+1;pnames{j} = sprintf('SunDisk2 Excess Pressure in Pa '); pscl(j)=1.0E-3;
    %     9 j=j+1;pnames{j} = sprintf('SunDisk2 Radius a in m         '); pscl(j)=1.0E3;
    %
    %     function [U1,U2,U3]=penny(xo,yo,xs,ys,nu,E,H,A,P)
    
    %kSunDisk = 46;
    kp = kSunDisk1 - koffset;
    if (abs(pg(kp+3)) > 0)
        switch idatatype1
            case 0 % observable is phase
                uSunDisk1 = penny(DST.x',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
            case -1 % observable is gradient
                tl = penny(DST.x'-DST.dx',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
                tr = penny(DST.x'+DST.dx',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
                uSunDisk1 = (tr - tl)/2.0;
            otherwise
                warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                uSunDisk1 = zeros(3,ndata);
        end
    else
        uSunDisk1 = zeros(3,ndata);
    end
    %kSunDisk = kSunDisk+5;
    kp = kSunDisk2 - koffset;
    if abs(pg(kp+3)) > 0
        switch idatatype1
            case 0 % observable is phase
                uSunDisk2 = penny(DST.x',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
            case -1 % observable is gradient
                tl = penny(DST.x'-DST.dx',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
                tr = penny(DST.x'+DST.dx',DST.y'...
                    ,pg(kp+0),pg(kp+1),nu,E ...
                    ,pg(kp+2),pg(kp+4),pg(kp+3));
                uSunDisk2 = (tr - tl)/2.0;
            otherwise
                warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
                uSunDisk2 = zeros(3,ndata);
        end
    else
        uSunDisk2 = zeros(3,ndata);
    end
    
    % 20121127
    %kcomsol = 67;
    kp = kCS - koffset;
%   pg(kCS+2)
    if (abs(pg(kp+2)) > 0)
        if numel(pt) ~= 2
            error('Cannot handle more than 2 epochs');
        end
        %for ipg=1:numel(pg)
        %    fprintf(1,'%4d %s %12.4e\n',ipg,char(PST.names{ipg+me*11}),pg(ipg));
        %end

        % parameter structure for Comsol
        %PSTCS = struct([]);
        secperyr = 365.25 * 24 * 3600;
        % initial epoch is 1980 - should be reference epoch eventually, but
        % using tquake triggers step function behavior
        PSTCS.t1 = (pt(1)-1980.0)*secperyr; % master epoch in seconds
        PSTCS.t2 = (pt(2)-1980.0)*secperyr; % slave epoch in seconds
        PSTCS.p1(1) = pt(1);
        PSTCS.p1(2) = pt(2);
        PSTCS.xcen = pg(kp+0); % Easting Coordinate
        PSTCS.ycen = pg(kp+1); % Northing Coordinate
%         PSTCS.POR  = pg(kp+2); % porosity
%         PSTCS.PER  = pg(kp+3); % permeability
        PSTCS.POR  = 10^pg(kp+2); % log10 porosity
        PSTCS.PER  = 10^pg(kp+3); % log10 permeability
        PSTCS.EYM = E;              % Young's modulus
        PSTCS.NUP = nu;             % Poisson's ratio
        PSTCS.verbose = 1;          % tell us what you are doing
        %PSTCS.verbose = 0;          % please be quiet
        %PSTCS.mname = '/data/iceland/Svart/T367_T410/comsol2d/svartsengi2DaxiV09'; % should be passed in datafile
        %PSTCS.mname = 'svartsengi2DaxiV09';
        %PSTCS.mname = '/data/bradys/TSX/T53/strip_008/comsol3d/bradys3dV03';
        PSTCS.mname = strrep(PST.datafilename,'.mph','');
        PSTCS.mparam = mparam;
        PSTCS.names = PST.names;
        % data structure for Comsol is the same
        %[times1,ucs,vcs,wcs] = funPoro2DaxiV09(PSTCS,DST);     % run from scratch
        [times1,xpts,ypts,zpts,ucs,vcs,wcs] = fun_comsol3d(PSTCS,DST); 
        % calculate range as dot product in meters
        UCS = -1.0* (ucs.*DST.uvx + vcs.*DST.uvy + wcs.*DST.uvz);
        
        clear ucs vcs wcs;
    else
        UCS = zeros(ndata,1);
    end
    
    % 20130521 adapt for abaqus
    %if abs(PST.p1(kAQ)) > 0.
    %         PSTAQ.mname = PST.datafilename;
    %         PSTAQ.pressure = PST.p1(kAQ+0); % Pressure increase in Pascals
    %         PSTAQ.xcen     = PST.p1(kAQ+1); % Easting coordinate of Center in meters
    %         PSTAQ.ycen     = PST.p1(kAQ+2); % Northing coordinate of Center in meters
    %         PSTAQ.depth    = PST.p1(kAQ+3); % Depth of center in meters
    %kAQ = kcomsol+4;
    kp = kAQ - koffset;
    if abs(pg(kp)) > 0.
        PSTAQ.mname    = PST.datafilename;
        PSTAQ.pressure = pg(kp+0); % Pressure increase in Pascals
%         PSTAQ.xcen     = pg(kp+1); % Easting coordinate of Center in meters
%         PSTAQ.ycen     = pg(kp+2); % Northing coordinate of Center in meters
%         PSTAQ.depth    = pg(kp+3); % Depth of center in meters
        PSTAQ.verbose  = 0;
        
%         % Run the model from scratch
%         switch idatatype1
%             case 0 % observable is phase
%                 uabaqus1a = abaqus2d_mogi1(PSTAQ,DST.x,DST.y);
%                 uabaqus1b = zeros(3,ndata);
%             case -1 % observable is gradient
%                 uabaqus1a = abaqus2d_mogi1(PSTAQ,DST.x-DST.dx,DST.y);
%                 uabaqus1b = abaqus2d_mogi1(PSTAQ,DST.x+DST.dx,DST.y);
%             otherwise
%                 warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
%                 uabaqus1a = zeros(3,ndata);
%                 uabaqus1b = zeros(3,ndata);
%         end
       
        % The model has already been run for a pressure of 1.0MPa
        scaleit = (PSTAQ.pressure / 1.0e6);
        uabaqus = (uabaqus1a - uabaqus1b) * scaleit;
    else
        uabaqus = zeros(3,ndata);
    end
    
    % radially dependent velocities
    %krv = kAQ+4;
    kp = kRad_Vel - koffset;
    if abs(pg(kp)) > 0. 
        %fprintf(1,'Starting Radial Velocity Model\n');
        % get time difference in years
        tdif_rv = DD * time_function(pt, 0.);       
        switch idatatype1
            case 0 % observable is phase
                % calculate radial velocity field
                Vrv = radial_velocity(pg(kp:kp+5),DST.x,DST.y);
                % convert velocities to displacements D Is R * T
                %Urv = Vrv * tdif_rv;
            case -1 % observable is gradient
                % calculate radial velocity field
                Vrvl = radial_velocity(pg(kp:kp+5),DST.x,DST.y);
                Vrvr = radial_velocity(pg(kp:kp+5),DST.x,DST.y);
                Vrv = Vrvr - Vrvl;
                % convert velocities to displacements D Is R * T
                %Urv = (Vrvr - Vrvl) * tdif_rv;
        end
        % position-dependent field of range change in meters
        %     rRad_Vel = -1.0* (colvec(Urv(1,:)).*DST.uvx ...
        %                 + colvec(Urv(2,:)).*DST.uvy ...
        %                 + colvec(Urv(3,:)).*DST.uvz);
        rRad_Vel = -1.0* tdif_rv .* (colvec(Vrv(1,:)).*DST.uvx ...
            +              colvec(Vrv(2,:)).*DST.uvy ...
            +              colvec(Vrv(3,:)).*DST.uvz);

    else
        %Urv = zeros(3,ndata);
        rRad_Vel = zeros(ndata,1);
    end
    
    % sum all RATES
    u = umogi1 + umogi2 + uokada1 + uokada2 + uokada3 + uYangPS + upinel + uSunDisk1 + uSunDisk2 + uabaqus;
    %disp 'u'; size u
    
    % position-dependent field of range change rate values in m/yr
    gmod = -1.0* (colvec(u(1,:)).*DST.uvx ...
                + colvec(u(2,:)).*DST.uvy ...
                + colvec(u(3,:)).*DST.uvz);
 
    %time function in a column vector
    %tdif = DD * colvec(pt);
    tdif = DD * time_function(pt, tquake);
    %fprintf(1,'Time difference in years %10.4f\n',tdif);
    
    % combine time and space dependence in a column vector
    rng1 = tdif .* gmod + nui + bas + UCS + rRad_Vel;  % unwrapped phase model in m

    % convert from meters to radians, returning a COLUMN vector
    rng = 2.0 * pi * rng1 ./ DST.mpercy;
    
    % range for deformation only, without accounting for nuisance parameters
    rng0 = tdif .* gmod * 2.0 * pi ./ DST.mpercy; % in radians
    
    %fprintf(1,'pair %d min = %12.4f max = %12.4f radians\n',i,min(min(rng1)),max(max(rng1)));
    
    
%     fprintf(1,'Minimal values in meters for xyzm %10.1f %10.1f %10.1f \n',min(DST.x),min(DST.y),min(DST.z));
%     fprintf(1,'Maximal values in meters for xyzm %10.1f %10.1f %10.1f \n',max(DST.x),max(DST.y),max(DST.z));
%     fprintf(1,'Minimal value  in meters for rRad_Vel %10.4f\n',min(rRad_Vel));
%     fprintf(1,'Maximal value  in meters for rRad_Vel %10.4f\n',max(rRad_Vel));
%     fprintf(1,'Minimal value  in meters for rng  %10.4f\n',min(rng));
%     fprintf(1,'Maximal value  in meters for rng  %10.4f\n',max(rng));
%     fprintf(1,'Minimal value  in meters for bas  %10.4f\n',min(bas));
%     fprintf(1,'Maximal value  in meters for bas  %10.4f\n',max(bas));
%     fprintf(1,'i    Time Function pt(i)\n');
%     for i = 1:numel(pt)
%         fprintf(1,'%5d %12.6f\n',i,pt(i));
%     end
    
    % return the modeled values
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

function ft = time_function(tepochs, tquake)
% return value of time function f(t)
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tquake  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch

[me, ncols] = size(tepochs);
if ncols == 1
    ft = zeros(me,1);
    
    if tquake >= nanmin(tepochs) && tquake <= nanmax(tepochs)
        %fprintf(1,'Using step function that turns on at epoch %f\n',tquake);
        itime=find(tepochs >= tquake);
        ft(itime) = 1.0;
    else
        ft = tepochs - tquake;
        %fprintf(1,'Linear in time (secular deformation) ft = %10.4f\n',ft);
    end
else
    ncols
    error('Dimension problem');
end
return



