function varargout = funfitbradythermal(varargin)
%function varargout = funfitbradythermal(varargin)


%% model for range change based on thermal contraction
% 20161101 Kurt Feigl

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

    % the following parameters are unchanged
    j=j+1;pnames{j} = sprintf('Poisson Ratio dimless          '); pscl(j)=1.0E-3;
    j=j+1;pnames{j} = sprintf('Shear Modulus in Pa            '); pscl(j)=1.0E10;
    j=j+1;pnames{j} = sprintf('Coeff of Thermal Expansion invK'); pscl(j)=1.0E-5;
    % 
    if numel(strfind(lower(char(datafilename)),'.mat')) > 0
        % data file exists
        if fexist(char(datafilename)) == 1
            j=j+1;pnames{j} = sprintf('TM Origin Easting in m          '); pscl(j)=1.0e3;
            j=j+1;pnames{j} = sprintf('TM Origin Northing in m         '); pscl(j)=1.0e3;
        end
    end
    

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
    
    for i = 1:ndata
        % Row corresponding to this pair
        DD( i,DST.kmast(i)) = -1; % subtract master
        DDM(i,DST.kmast(i)) = -1; % subtract master 20130629 - REPAIR
        DD( i,DST.kslav(i)) = +1; % from slave
        DDS(i,DST.kslav(i)) = +1; % from slave
    end
    idatatype1 = DST.idatatype(1);
    for i=1:ndata
        if DST.idatatype(i) ~= idatatype1
            fprintf(1,'idatatype differs for i = %d %d %d\n'...
                ,i,idatatype1,DST.idatatype(i));
        end
    end
    %incidence angle in radians from vertical
    bscale = acos(DST.uvz ./((DST.uvx).^2+(DST.uvy).^2+(DST.uvz).^2));
    bscale = bscale - mean(bscale);
    
    % do once to initialize
    upinel1a = zeros(3,ndata);
    upinel1b = zeros(3,ndata); 
    % Pointers for parameter names
    koffset   =  parameters_per_epoch * me;
  
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
        end
    end
    
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
     TST.parameters_per_epoch = parameters_per_epoch;
    TST.koffset    = koffset;
    TST.p0         = PST.p0;
    TST.p1         = PST.p1;
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
    koffset    = TST.koffset;
     
%     fprintf(1,'Unpacking parameter vector\n');

    %% get values of parameters
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
 
    for ikp = 1:numel(pg)
        fprintf(1,'%d %s %12.4e\n',ikp,PST.names{ikp+koffset},pg(ikp));
    end

    %% pointer index into array pg
    % 1 Poisson_Ratio_dimless___________   2.5000e-01
    % 2 Shear_Modulus_in_Pa_____________   3.0000e+10
    % 3 Coeff_of_Thermal_Expansion_invK_   0.0000e+00
    % 4 TM_Origin_Easting_in_m__________   0.0000e+00
    % 5 TM_Origin_Northing_in_m_________   0.0000e+00
    kp = 0;
    kp=kp+1; nu       = pg(kp);
    kp=kp+1; mu       = pg(kp);
    kp=kp+1; alphav   = pg(kp);
    kp=kp+1; tmx0     = pg(kp);
    kp=kp+1; tmy0     = pg(kp);
    
    %% check dimensions
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
    inan = find(isfinite(nui) == 0);
    if numel(inan) > 0
        warning(sprintf('zeroing %d nuisance values\n',numel(inan)));
        nui = zeros(size(DST.x0));
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
    
    %% read matlab file
    if numel(strfind(lower(char(PST.datafilename)),'.mat')) > 0
        % data file exists
        if fexist(char(PST.datafilename)) == 1
            S = load(char(PST.datafilename));
%             snames = fieldnames(S)
            %    Name           Size            Bytes  Class     Attributes
            %
            %   ncols          1x1                 8  double
            %   nrows          1x1                 8  double
            %   ux         12000x1             96000  double
            %   uy         12000x1             96000  double
            %   uz         12000x1             96000  double
            %   x          12000x1             96000  double
            %   y          12000x1             96000  double
            
            %% need times
            % itime1 =find(abs(timesteps_matlab-datetime(2011,12,24)) < days(100));
            % itime2 =find(abs(timesteps_matlab-datetime(2012,10,27)) < days(100));
            delta_t_mod = datetime(2012,10,27) - datetime(2011,12,24)
            dt_mod_yr = years(delta_t_mod);
            delta_t_obs = DST.tslav - DST.tmast;
            dt_obs_yr = years(delta_t_obs);
            
            
            %% interplolate onto query points where we have InSAR data
            UX = griddata(S.e,S.n,S.ue,DST.x,DST.y,'natural');
            UY = griddata(S.e,S.n,S.un,DST.x,DST.y,'natural');
            UZ = griddata(S.e,S.n,S.uv,DST.x,DST.y,'natural');
            
            % calculate range as dot product in meters and normalize to one
            % year
            rngtm = -1.0 * (UX.*DST.uvx + UY.*DST.uvy + UZ.*DST.uvz) .* dt_obs_yr / dt_mod_yr;
            
%             %% replace NaN with zeros!!
%             inan = find(isfinite(rngtm) == 0);
%             if numel(inan) > 0
%                 warning(sprintf('Replacing %d nan values with zeros\n',numel(inan)));
%                 rngtm(inan) = 0.;
%             end
%             
            %% replace outliers with zeros!!
            rstd = nanstd(rngtm);
            ravg = nanmean(rngtm);
            ibig = find(abs(rngtm-ravg) > 5*rstd);
            if numel(ibig) > 0
                warning(sprintf('Replacing %d big values with zeros\n',numel(ibig)));
                rngtm(ibig) = 0.;
            end

            
%             figure
%             histogram(rngtm);
%             title('rngtm');
%             clear UX UY UZ;
        else
            error(sprintf('Could not find data file named %s\n',PST.datafilename));
        end
    else
        rngtm = zeros(ndata,1);
    end
    
%     fprintf(1,'Summing vectors\n');

    % ************************************
    % DONE WITH INDIVIDUAL PARTS OF DISPLACEMENT FIELD
    % ************************************
%     % sum all RATES for separable part
%     u = umogi1 + umogi2 + uokada1 + uokada2 + uokada3 + uYangPS + upinel + uSunDisk1 + uSunDisk2;
% %     disp 'u'; size u
    
%     % position-dependent field of range change rate values in m/yr
%     gmod = -1.0 * (colvec(u(1,:)).*DST.uvx ...
%                 + colvec(u(2,:)).*DST.uvy ...
%                 + colvec(u(3,:)).*DST.uvz);
 
    %time function in a column vector
    %tdif = DD * colvec(pt);
    %tdif = DD * time_function(pt, tquake);
    %fprintf(1,'Time difference in years %10.4f\n',tdif);
%     tdif = DD * timefun(pt, tquake);
    %fprintf(1,'Time difference in years %10.4f\n',tdif);

    
    %% combine time and space dependence in a column vector
    %rng0 = tdif .* gmod + rngtm;  % modeled range change in m
    rng0 = rngtm;  % modeled range change in m
    
    %% add nuisance contribution
    
    rng1 = rng0 + nui + bas ;  

    %% create observable quantity 
    switch idatatype1
        case 0  % observable is phase or gradient in radians
            rng0 = 2.0 * pi * rng0 ./ DST.mpercy;
            rng =  2.0 * pi * rng1 ./ DST.mpercy;
       case -1  % observable is dimensionless gradient
            rng =  rng1;
        case 2  % observable is range change in meters
            rng = rng1;
        otherwise
            warning(sprintf('Unknown value of idatatype1 = %d\n',idatatype1));
            rng = zeros(ndata,1);
    end
    
%     % replace NaN with zeros?
%     inan = find(isfinite(rng0) == 0);
%     rng0(inan) = 0;
    
%    fprintf(1,'nanmin = %12.4f nanmax = %12.4f \n',nanmin(nanmin(rng1)),nanmax(nanmax(rng1)));
    
    
    fprintf(1,'Minimal values in meters for xyzm %10.1f %10.1f %10.1f \n',nanmin(DST.x),nanmin(DST.y),nanmin(DST.z));
    fprintf(1,'Maximal values in meters for xyzm %10.1f %10.1f %10.1f \n',nanmax(DST.x),nanmax(DST.y),nanmax(DST.z));
% %     fprintf(1,'Minimal value  in meters for rRad_Vel %10.4f\n',nanmin(rRad_Vel));
% %     fprintf(1,'Maximal value  in meters for rRad_Vel %10.4f\n',nanmax(rRad_Vel));
    fprintf(1,'Minimal value  in meters for rng  %10.4f\n',nanmin(rng));
    fprintf(1,'Maximal value  in meters for rng  %10.4f\n',nanmax(rng));
%     fprintf(1,'Minimal value  in meters for bas  %10.4f\n',nanmin(bas));
%     fprintf(1,'Maximal value  in meters for bas  %10.4f\n',nanmax(bas));
%     fprintf(1,'i    Time Function pt(i)\n');
%     for i = 1:numel(pt)
%         fprintf(1,'%5d %12.6f\n',i,pt(i));
%     end
    

%     fprintf(1,'Done calculating. Returning...\n');
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


