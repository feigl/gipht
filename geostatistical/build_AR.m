function AR = build_AR(A,d,data_dist,data_sigma,X,Y)
%% build A'*R where R is data spatial covariance matrix (too large to store)
% inputs 
% A == design matrix G in equation (12)
% d == data vector
% X == X coordinates of pixels in data vector [meters]
% Y == Y coordinates of pixels in data vector [meters]
% data_dist   == length scale for correlation function a = 230 m [equation
% (3)]
% data_sigma  == variance of data [m] \sigma_GPS = 0.005 m [equation (3)]


% Reference for equations
% Reinisch, E. C., M. Cardiff, and K. L. Feigl (2018), 
%Characterizing Volumetric Strain at Brady Hot Springs, Nevada, USA Using
%Geodetic Data, Numerical Models, and Prior Information,
%Geophysical Journal International, 1501?1513. 
%http://dx.doi.org/10.1093/gji/ggy347 
%The geothermal field at Brady Hot Springs, Nevada has subsided over the
%past decade. Between 2004 and 2014, the rate of downward vertical
%displacement was on the order of 10 millimeters per year, as measured by
%two independent geodetic techniques: Interferometric synthetic aperture
%radar (InSAR) and Global Positioning System (GPS). The observed
%deformation field forms an approximately elliptical bowl that is 4
%kilometers long and aligned with the trace of the NNE striking normal
%fault system. We use modeling to estimate the plausibility of pressure
%changes or thermal contraction as the cause of the observed subsidence. As
%a result, Bayesian inference favors with ?very strong evidence? thermal
%contraction over other hypotheses as the dominant driving mechanism for
%the observed subsidence. Using InSAR data spanning from 22 July 2016 and
%22 August 2017, we estimate the volume change rate in the significantly
%deforming volume to be (-29 +/- 3) thousand cubic meters per year and the
%total rate of change in thermal energy between -53 and -79 Megawatt. We
%infer the total volume of cubes where the estimated volumetric strain rate
%is significantly different from zero with 95 percent confidence to be 119
%million cubic meters. We find that the main region of significant cooling
%occurs between the injection and production well locations. This result
%supports the idea that highly permeable conduits along faults channel
%fluids from shallow aquifers to the deep reservoir tapped by the
%production wells.


% references:
% From: Kurt Feigl <feigl@wisc.edu>
% Date: Tuesday, December 15, 2020 at 2:44 PM
% To: "Reinisch, Elena Cristina" <ecreinisch@lanl.gov>
% Subject: Re: [EXTERNAL] Matlab code for Geostatistical inversion
% Would you please confirm that I am looking at the correct piece of code (attached)?
% https://raw.githubusercontent.com/feigl/gipht/Brady/HTC_mesh_build/mesh_range2d_HTC.m

% From: Reinisch, Elena Cristina <ecreinisch@lanl.gov>
% Sent: Tuesday, December 15, 2020 18:03
% To: Kurt Feigl <feigl@wisc.edu>
% Subject: Re: [EXTERNAL] Matlab code for Geostatistical inversion 
% That file is correct!
% values from original code:
%c = 0.01*(0.005/dt)^2; %3e-6; %5e-8;
%a_length_scale = 230; %m
% step_size = 3000

[arrows, arcols1] = size(A');
arcols = numel(d);
rowstep=1;

% handle scalar uncertainty
if numel(data_sigma) == 1
    data_sigma=data_sigma*ones(size(d));
else
    if numel(data_sigma) ~= numel(d)
        error('Length of uncertainty vector must match length of data vector\n');
end

if arcols1 ~= arcols
    arcols1
    arcols
    error(sprintf('Number of columns is incorrect.\n'));
else
    AR = zeros(size(A'));
    fprintf(1,'Starting %s...\n',mfilename);
    start_time = tic;
    
    % row index
    rind = [(0:rowstep:arcols), arcols];
    
    %% original code is a bit cryptic
    % rind = [(0:3000:arcols), arcols];
    % for i = 1:numel(rind)-1
    %     AR(:, rind(i)+1:rind(i+1)) = A'*[c*exp(-3.*(abs(sqrt((repmat(datx_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(datx_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2+(repmat(daty_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(daty_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2))./a))];
    % end
    %         ARfac = 1;
    %         AR = AR*ARfac;
    
    % this should be equivalent
    %whos
    nind = numel(rind)-1;
    for i = 1:nind
        X1 = repmat( X(rind(i)+1:rind(i+1))', [arcols, 1]);
        X2 = repmat( X(:)                   , [1, numel(rind(i)+1:rind(i+1))]);
        Y1 = repmat( Y(rind(i)+1:rind(i+1))', [arcols, 1]);
        Y2 = repmat( Y(:)                   , [1, numel(rind(i)+1:rind(i+1))]);
        T1 = exp( -3.*(abs(sqrt( (X1 - X2).^2 + (Y1 - Y2).^2)/data_dist))); % Equation (3)
        %P1 = (data_sigma^2)*(A'*T1);  % product
        % 2021/09/07 Kurt
        j=min([i,numel(data_sigma)]);
        P1 = (data_sigma(j)^2) * (A'*T1);  % product
        %     whos
        %     disp('c');size(c)
        %     disp('At');size(A')
        %     disp('T1');size(T1)
        %     disp('P1');size(P1)
        %     disp('size(AR(:, rind(i)+1:rind(i+1)))');size(AR(:, rind(i)+1:rind(i+1)))
        AR(:, rind(i)+1:rind(i+1)) = P1;
%         if mod(i,10) == 1
%             fprintf(1,'i = %6d of %6d Elapsed time in seconds: %.1f\n',i,nind,toc(start_time));
%         end
    end
end


return
end

