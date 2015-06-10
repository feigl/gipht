function varargout = fun_comsol3d(varargin)
%function [times,xpts,ypts,zpts,dupts,dvpts,dwpts] = fun_comsol3d(PST,DST,TST)
%varargout = funseparable27(varargin)
% number of arguments entered as input

% number of arguments entered as input
nin = nargin

% number of arguments to return as output
nout = nargout

% decide what to do based on number of arguments

if nin >= 1 PST = varargin{1}; end
if nin >= 2 DST = varargin{2}; end
if nin >= 3 TST = varargin{3}; end


savefile = sprintf('%s.mat',PST.mname);

clockstr = clock;
runname=sprintf('%04d%02d%02d_%02d%02d%02d'...
    ,clockstr(1),clockstr(2),clockstr(3),clockstr(4),clockstr(5),round(clockstr(6)));


% number of epochs
me = 0;
for i=1:PST.mparam
    if strfind(PST.names{i},'time_fn') > 0
        me = me+1;
    end
end


if me ~= 2
    error('Can only handle 2 epochs now.');
end

secperyr = 3600*24*365.25 ; % seconds per year

% just get the names of the parameters
if ismember(nin, [1 2]) == 1
    if PST.verbose == 1;
        fprintf(1,'Starting %s\n',PST.mname);
    end;
    
    %model = svartsengi2DaxiV09a;
    
    model = mphload(sprintf('%s.mph',PST.mname));
    
    param_names = model.param.varnames;
    for ip = 1:numel(param_names)
        pcode = char(param_names(ip));
        descr = char(model.param.descr(param_names(ip)));
        value = char(model.param.get(param_names(ip)));
        
        if PST.verbose == 1;
            fprintf(1,'    %-8s %-32s %s\n',pcode,value,descr);
        end
        
        % match parameter names here
        if findstr(pcode,'XMIN') > 0
            %XMIN = str2num(value)
            XMIN = sscanf(value,'%f')
        end
        if findstr(pcode,'YMIN') > 0
            %YMIN = str2num(value)
            YMIN = sscanf(value,'%f')
        end
    end
end

% run the model if needed
% This looks anywhere in the Matlab command search path
% if exist(savefile,'file')  ~= 2
% TST is empty
%if exist('TST','var') == 0

if ismember(nin, [2]) == 1
    if size(DST.x) ~= size(DST.y)
        error('size problem');
    end
    
    %     %     if PST.verbose == 1; mphmodel(model); end;
    
    %     % Defaults as of 2012NOV26 to produce about 40 mm of subsidence
    %     model.param.set('DTS', '30[year]', 'Delta Time Seconds');
    %     model.param.set('VIS', '2.84e-4[Pa*s]', 'Dynamic Viscosity of Fluid');
    %     %    model.param.set('POR', '0.2', 'Porosity');
    %     %   model.param.set('PER', '1.e-15[m^2]', 'Permeability');
    %     model.param.set('QMF', '-220[kg/s]', 'Q Mass Flux');
    %     model.param.set('CMF', '1/2.2e9[1/Pa]', 'Compressibility of Fluid');
    %     model.param.set('CMS', '4e-7[1/Pa]', 'Compressibility of Steam');
    %     model.param.set('HH0', '0[m]', 'Hydraulic Head at t=0 and z=0');
    %     %   model.param.set('EYM', '5[GPa]', 'E Young''s Modulus');
    %     %   model.param.set('NUP', '0.27', 'Poisson''s Ratio');
    %     model.param.set('DBC', '3*(1-2*NUP)/EYM [1/Pa]', 'Drained Bulk Compressibility');
    %     model.param.set('BOR', '10[m]', 'Radius of Borehole');
    %     model.param.set('RAD', '50e3[m]', 'Radius of Problem Domain');
    %     model.param.set('CAP', '400 [m]', 'Cap Rock Thickness');
    %     model.param.set('SFZ', '300 [m]', 'Shallow Feader Zone Thickness');
    %     model.param.set('BOD', '5000 [m]', 'Bottom depth');
    %     model.param.set('FZT', '800[m]', 'Feeder Zone Thickness');
    %
    %     % set Poisson's ratio
    %     model.param.set('NUP', sprintf('%e',PST.NUP),'Poisson''s Ratio');
    %
    %     % set Porosity
    %     model.param.set('POR', sprintf('%e',PST.POR),'Porosity in aquifer');
    %
    %     % set Permeability
    %     model.param.set('PER', sprintf('%e[m^2]',PST.PER),'Permeability');
    %
    %     % set Young's Modulus
    %     model.param.set('EYM', sprintf('%e[Pa]',PST.EYM), 'Young''s Modulus');
    
    % Try top decide if we have to run solver
    %prop=mphgetproperties(model.result);
    %if numel(prop) == 0
    try
        info=mphsolinfo(model,'NU','on')
    catch
        %if info.NUsol < 1
        %       % run the Darcy's Law part of the problem
        %     % Run the Poroelastic part of the problem
        %     % do not need the next line if the .mph file already contains nodeput
        if PST.verbose == 1; fprintf(1,'Starting to calculate %s\n',PST.mname); end;
        tic1=tic;
        model.sol('sol1').runAll;
        %     %model.sol('sol2').runAll; % not needed if Darcy and Poro solutions are combined
        %    model.sol('sol1').updateSolution; % no need to remesh
        
        mphsave(model,sprintf('%s.mph',PST.mname));
        if PST.verbose == 1; fprintf(1,'Finished calculating %s in %.f seconds\n',PST.mname,toc(tic1));end
    end
    
    info = mphsolinfo(model);
    if PST.verbose == 1
        info
    end
    
    
    % Comsol evaluation times in seconds
    times = info.solvals;
end

% now evaluate the output
if ismember(nin, [2 3]) == 1
    % start (master) and stop (slave) epochs for interferometric pairs
    % in seconds
    tsec1 =  0.;
    tsec2 = (PST.p1(2)- PST.p1(1))* secperyr;
    % FIX: should be time elapsed from some reference time.
    
    % obtain nodal coordinates first:
    nodestruct = mphxmeshinfo(model,'soltag','sol1');
    coords = nodestruct.nodes.coords;
    [ndim,ncols] = size(coords);
    x_node = coords(1,:); % x coordinate  in meters, w.r.t. xcen
    y_node = coords(2,:); % y coordinate  in meters, w.r.t. ycen
    z_node = coords(3,:); % z component in meters, positive upwards, w.r.t. 0
    
    x_node = x_node + XMIN;
    y_node = y_node + YMIN;
    %%% z_node = z_node + ZMIN;
    
    if PST.verbose == 1;
        figure;
        plot(x_node,y_node,'k+');
        xlabel('X-coordinate');
        ylabel('Y-coordinate');
        title('x_node,y_node');
        
        figure;
        %         plot(x_node,z_node,'b+');
        xlabel('X-coordinate');
        ylabel('Z-coordinate');
        title('x_node,z_node');
    end
    
    
    % pressure
    %p_node = mphinterp(model,{'p'},'coord',coords,'t',[PST.t1 PST.t2]); % pressure in Pa
    
    % displacement
    u_node = mphinterp(model,{'u'},'coord',coords,'t',[tsec1 tsec2]); % easting component
    v_node = mphinterp(model,{'v'},'coord',coords,'t',[tsec1 tsec2]); % northing component
    w_node = mphinterp(model,{'w'},'coord',coords,'t',[tsec1 tsec2]); % vertical component
    
    % stresses at nodes
    %Sxy_node = mphinterp(model,{'poro.sxy'},'coord',coords,'t',[tsec1 tsec2]); % xy component
    
    if PST.verbose == 1;
        fprintf(1,'Saving %s\n',savefile);
    end;
    %save(savefile, 'mname','runname','HHR','POR','times','r_node','z_node','u_node','w_node');
    %save(savefile,'runname','PST','times','r_node','z_node','u_node','w_node');
    clear model nodestruct coords info iset1 iset2;
    
    % find indices to times for evaluation
    [ntimes,n_nodes] = size(w_node);
    % if PST.verbose == 1
    %     ntimes
    %     n_nodes
    % end
    
    % it1 = find(abs(times-PST.t1)/secperyr < 0.5); % t1 = start time at 1 year
    % it2 = find(abs(times-PST.t2)/secperyr < 0.5); % t2 = stop  time at 2 years
    % it1 = unique(min(it1));
    % it2 = unique(max(it2));
    it1 = 1;
    it2 = 2;
    
    % get indices of coordinates on top
    itop = find(abs(z_node) < 1.0);
    
    %itop = find(z_node >= max(max(DST.z)));
    
    if numel(itop) < 10
        itop
        error('Not enough nodes on top to interpolate');
    end
    % x_node = colvec(x_node(itop));
    % y_node = colvec(y_node(itop));
    % z_node = colvec(z_node(itop));
    % u_node = u_node(:,itop);
    % v_node = v_node(:,itop);
    % w_node = w_node(:,itop);
    
    
    
    % rate  during time interval from t1 to t2
    % dhdt = secperyr*(u_node(it2,:) - u_node(it1,:))/(times(it2)-times(it1));
    % dwdt = secperyr*(w_node(it2,:) - w_node(it1,:))/(times(it2)-times(it1));
    
    % incremental displacement during time interval from t1 to t2
    du_node = colvec(u_node(it2,:) - u_node(it1,:));
    dv_node = colvec(v_node(it2,:) - v_node(it1,:));
    dw_node = colvec(w_node(it2,:) - w_node(it1,:));
    
    
    if PST.verbose == 1;
        figure
        hist(colvec(du_node));
        xlabel('du_node');
        ylabel('N');
    end
    
    % incremental change in stress during time interval from t1 to t2
    %dSxy_node = colvec(Sxy_node(it2,:) - Sxy_node(it1,:));
    
    % evaluate functions at observation locations
    % Vq = interp2(X,Y,V,Xq,Yq) interpolates to find Vq, the values of the
    %     underlying 2-D function V at the query points in matrices Xq and Yq.
    %     Matrices X and Y specify the points at which the data V is given.
    %  dupts  = interp2(x_node,y_node,du_node,DST.x,DST.y,'cubic'); % easting component
    %  dvpts  = interp2(x_node,y_node,dv_node,DST.x,DST.y,'cubic'); % northing component
    %  dwpts  = interp2(x_node,y_node,dw_node,DST.x,DST.y,'cubic'); % vertical component
    
    % evaluate functions at observation locations
    % create interpolating functions
    % Fu = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(du_node(itop)),'natural'); % easting component
    % Fv = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dv_node(itop)),'natural'); % northing component
    % Fw = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dw_node(itop)),'natural'); % vertical component
    
    
    % % plot values at _nodes of mesh
    % xpts = x_node;
    % ypts = y_node;
    % zpts = z_node;
    
    % plot values at observation points
    xpts = DST.x;
    ypts = DST.y;
    zpts = DST.z;
    
    % evaluate interpolating functions at requested locations
    % dupts  = Fu(xpts,ypts); % horizontal component
    % dvpts  = Fv(xpts,ypts); % horizontal component
    % dwpts  = Fw(xpts,ypts); % horizontal component
    
    % estimate interpolant functions and evaluate them
    dupts = griddata2(colvec(x_node(itop)),colvec(y_node(itop)),colvec(du_node(itop)),xpts,ypts);
    dvpts = griddata2(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dv_node(itop)),xpts,ypts);
    dwpts = griddata2(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dw_node(itop)),xpts,ypts);
    
    
    if PST.verbose == 1;
        figure
        hist(colvec(dupts));
        xlabel('dupts');
        ylabel('N');
        
        figure
        whos xpts ypts dupts dvpts
        quiver(xpts,ypts,dupts,dvpts);
        title(sprintf('displacements from %s',PST.mname));
        xlabel('X coordinate');
        ylabel('Y coordinate');
        
        % make profile
        figure; hold on;
        yprof = nanmean(y_node);
        iprof = find(abs(y_node(itop) - yprof) < 1000.);
        
        
        if max(iprof) > numel(dupts)
            iprof
            whos iprof
            whos dupts
            error('This should not happen.');
        end
        plot(colvec(x_node(iprof)),colvec(dupts(iprof)),'ro');
        plot(colvec(x_node(iprof)),colvec(dvpts(iprof)),'b^');
        plot(colvec(x_node(iprof)),colvec(dwpts(iprof)),'ks');
        legend('dupts','dvpts','dwpts');
        title(sprintf('Displacements from %s',PST.mname));
        xlabel('X coordinate');
        ylabel('Displacement [m]');
        
        %     figure;
        %     yprof = nanmean(y_node);
        %     iprof = find(abs(y_node - yprof) < 10.);
        %     plot(x_node(iprof),dSxy_node(iprof)/1.e6,'ro');
        %     title(sprintf('stress from %s',PST.mname));
        %     xlabel('X coordinate');
        %     ylabel('stress (MPa)');
        %
    end
    % else
    %     %     if PST.verbose == 1; fprintf(1,'Loading %s\n',savefile);end
    %     %     load(savefile);
    
    % return output variables to calling routine for further use
    % times,xpts,ypts,zpts,dupts,dvpts,dwpts
    varargout(1) = {times};
    varargout(2) = {xpts};
    varargout(3) = {ypts};
    varargout(4) = {zpts};
    varargout(5) = {dupts};
    varargout(6) = {dvpts};
    varargout(7) = {dwpts};    
end

return
end

