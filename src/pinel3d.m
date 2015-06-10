function varargout = pinel3d(varargin)
%function u = pinel3c(xinm,yinm,g,E,nu,rho,loadfile)
% calculate displacements [east, north, up] at [xobs, yobs] due to a of height HI on grid XI, YI
% Pinel et al. (2008) Geophys. J. Int. v. 169, pp 325-338
% assume value of density (rho) is same everywhere
% Kurt Feigl 2010-OCT-05
% number of arguments entered as input



nin = nargin;

% number of arguments to return as output
nout = nargout;



% decide what to do based on number of arguments
if nin == 1 && nout == 4
    % [XI,YI,HI,mn] = pinel3d(char(PST.datafilename));
    % just make grids
    
    % loadfile name
    datafilename = varargin{1};
    
    if fexist(datafilename) > 0
        fprintf(1,'Starting %s with datafilename = %s\n'...
            ,mfilename,datafilename);
    else
        error(sprintf('Could not find datafilename = %s\n'...
            ,datafilename));
    end
    
    matfilename = sprintf('%s.mat',datafilename);
    if fexist(matfilename) > 0
        load(matfilename,'-mat');
        if strcmp(datafilenamesave,datafilename) == 1
                    fprintf(1,'In %s using existing grids from matfilename = %s\n'...
            ,mfilename,matfilename);
        else
            error('Data file names do not match');
        end
    else
         fprintf(1,'In %s making new grids from datafilename = %s\n'...
            ,mfilename,datafilename);
        %loadfile = 'VaWEmism93-95.dat'
        L = load(datafilename);
        xl = L(:,1);
        yl = L(:,2);
        hl = L(:,3);
        
        % grid spacing in meters
        if numel(strfind(lower(datafilename),'lava')) > 0
            % grid spacing in meters - for Okmok lava
            dX = 50;
            dY = 50;
        else
            dX = 1000;
            dY = 1000;
        end
        
        xax = [min(xl):dX:max(xl)];
        yax = [min(yl):dY:max(yl)];
        
        [XI, YI] = meshgrid(xax,yax);
        %HI = griddata(xl,yl,hl,XI,YI,'cubic'); % makes NaN
        %HI = griddata(xl,yl,hl,XI,YI,'v4'); % unreliable
        HI = griddata2(xl,yl,hl,XI,YI,'nearnat'); % 20130701
        
        mn = size(HI);
        
        % fill holes with zeros
        HI(isnan(HI))= 0;
        % fill negative values with zero for Okmok case
        if numel(strfind(lower(datafilename),'lava')) > 0
            HI(find(HI<=0.0)) = 0;
        end
        
        %disp 'number of nonzero values'
        %numel(find(abs(HI) > 0))
        %     figure;hold on;
        %     %imagesc(xax/1000,yax/1000,HI);
        %     imagesc(HI)
        %     xlabel('Easting (km)');
        %     ylabel('Northing (km)');
        %     axis xy; axis tight; axis equal;
        %     colorbar;
        %     title(datafilename);
        %     %feval(printfun,sprintf('%s.pdf',datafilename));
        
        datafilenamesave = datafilename;
        save(matfilename,'XI','YI','HI','mn','datafilenamesave');
    end
    
    
    % pass variables back to calling routine for future use
    varargout(1) = {XI};
    varargout(2) = {YI};
    varargout(3) = {HI};
    varargout(4) = {mn};
    return
end


if nin == 6 && nout == 1
    % calculate unit dispalcements
    %u = pinel3d(DST.x,DST.y,XI,YI,HI,mn);
    %              1     2   3, 4,  5,6)
    xobs    = varargin{1};
    yobs    = varargin{2};
    XI      = varargin{3};
    YI      = varargin{4};
    HI      = varargin{5};
    mn      = varargin{6};
    
    %     XI=reshape(XI,mn(1),mn(2));
    %     YI=reshape(YI,mn(1),mn(2));
    %     HI=reshape(HI,mn(1),mn(2));
    
    %     disp XI;size(XI)
    %     disp YI;size(YI)
    %     disp HI;size(HI)
    
    % constant part of eq (1)
    %ur = -1*rho*g*(1+nu)*(1-2*nu)/(2*pi*E);
    ur = 1.0;
    
    % constant part of eq (2)
    %uz = g * rho* (1- nu^2) / (pi * E);
    uz = 1.0;
    
    % initialize
    ndata = numel(xobs);
    u = zeros(3,ndata);
    
    % increments
    dx = XI(1,2)-XI(1,1);
    dy = YI(2,1)-YI(1,1);
    
    % calculate displacement at each observation point as integral over all grid points
    for i=1:ndata
        % position of obs point w.r.t. grid point
        x = xobs(i) - XI;
        y = yobs(i) - YI;
        r = sqrt(x.^2 + y.^2);
        th = angle(complex(x./r,y./r));
        % displacements
        u(1,i) =      sum(sum(cos(th).*HI .* ur .* dx .* dy./r));
        u(2,i) =      sum(sum(sin(th).*HI .* ur .* dx .* dy./r));
        u(3,i) = -1.0*sum(sum(         HI .* uz .* dx .* dy./r));
    end
    
    % trim out NaN
    inan = find(isfinite(u)==0);
    u(inan) = 0;
    % return the displacements
    varargout(1) = {u};
    return
end

% just scale unit displacements by new values of material properties
if nin == 5 && nout == 1
    
    %u = pinel3d(g,E,nu,rho,u1);
    %            1,2,3 ,4,  5
    g       = varargin{1};
    E       = varargin{2};
    nu      = varargin{3};
    rho     = varargin{4};
    % velocity field calculated with ur = uz = 1;
    u1      = varargin{5};
    %u1      = reshape(u1,3,numel(u1)/3);
    %     disp u1;size(u1)
    %     max(max(u1(1,:)))
    %     max(max(u1(2,:)))
    %     max(max(u1(3,:)))
    
    % constant part of eq (1)
    ur = -1*rho*g*(1+nu)*(1-2*nu)/(2*pi*E);
    
    % constant part of eq (2)
    uz = g * rho* (1- nu^2) / (pi * E);
    
    % scale displacments
    u(1,:) = ur * u1(1,:);
    u(2,:) = ur * u1(2,:);
    u(3,:) = uz * u1(3,:);
    
    % return the displacements
    varargout(1) = {u};
    return
end





