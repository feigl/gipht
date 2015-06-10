function [dupts,dvpts,dwpts] = fun_abaqus2d(PST,DST,TST)
% calculate displacements using Abaqus
% 2013-05-21 Kurt Feigl

%mname = 'svartsengi2DaxiV01';
%mname = PST.mname;

if numel(DST.x) ~= numel(DST.y)
    error('size problem');
end

savefile = sprintf('%s.mat',PST.mname);

clockstr = clock;
runname=sprintf('%04d%02d%02d_%02d%02d%02d'...
    ,clockstr(1),clockstr(2),clockstr(3),clockstr(4),clockstr(5),round(clockstr(6)));

% number of arguments entered as input
nin = nargin;

%secperyr = 3600*24*365.25 ; % seconds per year

% run the model if needed
if nin == 2
    if PST.verbose == 1;
        fprintf(1,'Starting %s\n',PST.mname);
    end;

%     if PST.verbose == 1
%         info
%     end
%     
%     
%     if PST.verbose == 1;
%         param_names = model.param.varnames;
%         for ip = 1:numel(param_names)
%             pcode = char(param_names(ip));
%             descr = char(model.param.descr(param_names(ip)));
%             value = char(model.param.get(param_names(ip)));
%             
%             % could match parameter names here
%             fprintf(1,'    %-8s %-32s %s\n',pcode,value,descr);
%         end
%     end
%     
%     % evaluation times in seconds
%     times = info.solvals;
%     
%     % obtain nodal coordinates first:
%     nodestruct = mphxmeshinfo(model,'soltag','sol1');
%     coords = nodestruct.nodes.coords;
%     [ndim,ncols] = size(coords);
%     x_node = coords(1,:); % x coordinate  in meters, w.r.t. xcen
%     y_node = coords(2,:); % y coordinate  in meters, w.r.t. ycen
%     z_node = coords(3,:); % z component in meters, positive upwards, w.r.t. 0
%     
%     % pressure
%     %p_node = mphinterp(model,{'p'},'coord',coords,'t',[PST.t1 PST.t2]); % pressure in Pa
%     
%     % displacement
%     u_node = mphinterp(model,{'u'},'coord',coords,'t',[PST.t1 PST.t2]); % easting component
%     v_node = mphinterp(model,{'v'},'coord',coords,'t',[PST.t1 PST.t2]); % northing component
%     w_node = mphinterp(model,{'w'},'coord',coords,'t',[PST.t1 PST.t2]); % vertical component
%     
%     % stresses at nodes
%     Sxy_node = mphinterp(model,{'poro.sxy'},'coord',coords,'t',[PST.t1 PST.t2]); % xy component


    convert_abaqus_dat_to_csv('mogi-test.dat');
    [tstep,rnode,znode,unode,wnode] = read_abaqus_csv2d('mogi-test.csv');
    
    % interpolate the 2-D field to 3-D by cylindrical symmetry
    [upts,vpts,wpts] = interpolate2Ddisp(PST,DST,rnode,znode,unode,wnode);
    mpts = sqrt(upts.^2 + vpts.^2 + wpts.^2);


    if PST.verbose == 1;
        fprintf(1,'Saving %s\n',savefile);
    end;
    %save(savefile, 'mname','runname','HHR','POR','times','r_node','z_node','u_node','w_node');
    %save(savefile,'runname','PST','times','r_node','z_node','u_node','w_node');
    clear model nodestruct coords info iset1 iset2;
else
    %     if PST.verbose == 1; fprintf(1,'Loading %s\n',savefile);end
    %     load(savefile);
end



% find indices to times for evaluation
%[ntimes,n_nodes] = size(w_node);
% if PST.verbose == 1
%     ntimes
%     n_nodes
% end

% it1 = find(abs(times-PST.t1)/secperyr < 0.5); % t1 = start time at 1 year
% it2 = find(abs(times-PST.t2)/secperyr < 0.5); % t2 = stop  time at 2 years
% it1 = unique(min(it1));
% it2 = unique(max(it2));
%it1 = 1;
%it2 = 2;

% get indices of coordinates on top
%itop = find(abs(z_node) < 1.0);
%itop = find(z_node >= max(max(DST.z)));
% x_node = colvec(x_node(itop));
% y_node = colvec(y_node(itop));
% z_node = colvec(z_node(itop));
% u_node = u_node(:,itop);
% v_node = v_node(:,itop);
% w_node = w_node(:,itop);



if PST.verbose == 1;
    % plot the results
    figure; hold on;
    quiver(DST.x/1.e3,DST.y/1.e3,upts,vpts);
    xlabel('Easting coordinate [km]');
    ylabel('Northing coordinate [km]');
    
    figure;hold on;
    iok = find(abs(DST.y) <= 1.);
    plot(DST.x(iok)/1.e3,upts(iok)*1.e3,'r+-');
    plot(DST.x(iok)/1.e3,vpts(iok)*1.e3,'go-');
    plot(DST.x(iok)/1.e3,wpts(iok)*1.e3,'b^-');
    plot(DST.x(iok)/1.e3,mpts(iok)*1.e3,'k*-');
    xlabel('Easting coordinate [km]');
    ylabel('Displacement [mm]');
    legend('Eastward','Northward','Upward','Magnitude');
    title('mogi-test');
    
end

% rate  during time interval from t1 to t2
% dhdt = secperyr*(u_node(it2,:) - u_node(it1,:))/(times(it2)-times(it1));
% dwdt = secperyr*(w_node(it2,:) - w_node(it1,:))/(times(it2)-times(it1));

% % incremental displacement during time interval from t1 to t2
% % du_node = colvec(u_node(it2,:) - u_node(it1,:));
% % dv_node = colvec(v_node(it2,:) - v_node(it1,:));
% % dw_node = colvec(w_node(it2,:) - w_node(it1,:));
% % 
% % % incremental change in stress during time interval from t1 to t2
% % dSxy_node = colvec(Sxy_node(it2,:) - Sxy_node(it1,:));
% 
% % evaluate functions at observation locations
% % Vq = interp2(X,Y,V,Xq,Yq) interpolates to find Vq, the values of the
% %     underlying 2-D function V at the query points in matrices Xq and Yq.
% %     Matrices X and Y specify the points at which the data V is given.
% %  dupts  = interp2(x_node,y_node,du_node,DST.x,DST.y,'cubic'); % easting component
% %  dvpts  = interp2(x_node,y_node,dv_node,DST.x,DST.y,'cubic'); % northing component
% %  dwpts  = interp2(x_node,y_node,dw_node,DST.x,DST.y,'cubic'); % vertical component
% 
% % evaluate functions at observation locations
% % create interpolating functions
% % Fu = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(du_node(itop)),'natural'); % easting component
% % Fv = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dv_node(itop)),'natural'); % northing component
% % Fw = TriScatteredInterp(colvec(x_node(itop)),colvec(y_node(itop)),colvec(dw_node(itop)),'natural'); % vertical component
% 
% 
% % % plot values at _nodes of mesh
% % xpts = x_node;
% % ypts = y_node;
% % zpts = z_node;
% 
% plot values at observation points
xpts = DST.x;
ypts = DST.y;
zpts = DST.z;
% 
% % evaluate interpolating functions at requested locations
% dupts  = Fu(xpts,ypts); % horizontal component
% dvpts  = Fv(xpts,ypts); % horizontal component
% dwpts  = Fw(xpts,ypts); % horizontal component
% no time dependence yet
dupts = upts;
dvpts = vpts;
dwpts = wpts;

% 
% if PST.verbose == 1;
%     figure
%     quiver(xpts,ypts,dupts,dvpts);
%     title(sprintf('displacements from %s',PST.mname));
%     xlabel('X coordinate');
%     ylabel('Y coordinate');
%     
%     figure;
%     iprof = find(abs(y_node - 4406.5e3) < 1.);
%     plot(x_node(iprof),dSxy_node(iprof)/1.e6,'ro');
%     title(sprintf('stress from %s',PST.mname));
%     xlabel('X coordinate');
%     ylabel('stress (MPa)');
%     
% end
return

