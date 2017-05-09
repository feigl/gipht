function plot_temperature_field(PTFIELD,GRID,ORIGIN,p_or_t,dodiffs,date_times)
%% plot temperature field
% 20170131 Kurt Feigl

% initialize
nf = 0;
[nvoxelsG,dummy] = size(GRID.GridBlock);

%% Find indices of time steps
ktimes = zeros(numel(date_times),1);
for i = 1:numel(date_times)
    ndays = abs(days((PTFIELD.t - date_times(i))))
    ndays_min = unique(nanmin(ndays))
    kt=find(abs(ndays - ndays_min) <= eps);
    ktimes(i) = kt(1)
end

%% time interval
if dodiffs == 1 && numel(date_times) == 2
    dt = days(PTFIELD.t(ktimes(2)) - PTFIELD.t(ktimes(1)));
else
    dt = 0;
end


%% find some statistics
Tmin = min(colvec(PTFIELD.T))
Tmax = max(colvec(PTFIELD.T))
Tmed = median(colvec(PTFIELD.T))
Tstd = std(colvec(PTFIELD.T))

Pmin = min(colvec(PTFIELD.P))
Pmax = max(colvec(PTFIELD.P))
Pmed = median(colvec(PTFIELD.P))
Pstd = std(colvec(PTFIELD.P))


%% find corresponding indices
[iTmin,kTmin] = find(abs(PTFIELD.T - Tmin) <= eps);iTmin=iTmin(1);
[iTmax,kTmax] = find(abs(PTFIELD.T - Tmax) <= eps);iTmax=iTmax(1);
%[iTmed,kTmed] = find(abs(PTFIELD.T - Tmed) <= eps);iTmed=iTmed(1);
[dummy,iTmed] = sort(abs(PTFIELD.T - Tmed));iTmed=iTmed(1);

iTmaxstd=find(abs(std(PTFIELD.T')-max(std(PTFIELD.T'))) <= eps);iTmaxstd=iTmaxstd(1);

[iPmin,kPmin] = find(abs(PTFIELD.P - Pmin) <= eps);iPmin=iPmin(1);
[iPmax,kPmax] = find(abs(PTFIELD.P - Pmax) <= eps);iPmax=iPmax(1);
% [iPmed,kPmed] = find(abs(PTFIELD.P - Pmed) <= eps);iPmed=iPmed(1);
[dummy,iPmed] = sort(abs(PTFIELD.P - Pmed));iPmed=iPmed(1);

iPmaxstd=find(abs(std(PTFIELD.P')-max(std(PTFIELD.P'))) <= eps);iPmaxstd=iPmaxstd(1);

tzero = datetime(0,1,1);
tzero.TimeZone = 'UTC';




%% choose location for slices
% xcen = mean([min(GRID.ModelX) max(GRID.ModelX)])
% ycen = mean([min(GRID.ModelY) max(GRID.ModelY)])
% zcen = mean([min(GRID.ModelZ) max(GRID.ModelZ)])

xcen=GRID.ModelX(iTmax)
ycen=GRID.ModelY(iTmax)
zcen=GRID.ModelZ(iTmax)
ecen=GRID.Easting(iTmax)
ncen=GRID.Northing(iTmax)
vcen=GRID.Elevation(iTmax)
%
%% plot time series of Temperature for various voxels
nf=nf+1;h(nf)=figure;hold on;
plot(years(PTFIELD.t-tzero),PTFIELD.T(iTmin,1:end),'b-');
plot(years(PTFIELD.t-tzero),PTFIELD.T(iTmax,1:end),'r-');
plot(years(PTFIELD.t-tzero),PTFIELD.T(iTmed,1:end),'k-');
plot(years(PTFIELD.t-tzero),PTFIELD.T(iTmaxstd,1:end),'g-');
legend('min','max','median','max std dev');
xlabel('time step [year]');
ylabel('Temperature [deg C]');
title(PTFIELD.dirname,'Interpreter','none');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));


%% plot time series of Pressure for various voxels
nf=nf+1;h(nf)=figure;hold on;
plot(years(PTFIELD.t-tzero),PTFIELD.P(iPmin,1:end),'b-');
plot(years(PTFIELD.t-tzero),PTFIELD.P(iPmax,1:end),'r-');
plot(years(PTFIELD.t-tzero),PTFIELD.P(iPmed,1:end),'k-');
plot(years(PTFIELD.t-tzero),PTFIELD.P(iPmaxstd,1:end),'g-');
%xlabel(sprintf('time step [days after %s]',strt0));
legend('min','max','median','max std dev');
xlabel('time step [year]');
ylabel('Pressure [Pa]');
title(PTFIELD.dirname,'Interpreter','none');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));


%% draw voxels
% faces are defined by vertex number
%https://ths1104geek.wordpress.com/2010/08/25/draw-a-polyhedron-with-matlab/
for j = 1:3
    % draw voxels
    
    % select voxels
    switch j
        case 1 % cross-section in vertical plane normal to X-axis
            iplot = intersect(1:nvoxelsG,find(abs(GRID.ModelX - xcen) <= GRID.ModelDX_m_/2.));
            korners = [1,4,8,5];
        case 2 % cross-section in vertical plane normal to Y-axis
            iplot = intersect(1:nvoxelsG,find(abs(GRID.ModelY - ycen) <= GRID.ModelDY_m_/2.));
            korners = [4,3,7,8];
        case 3 % map view in horizontal plane normal to Z-axis
            iplot = intersect(1:nvoxelsG,find(abs(GRID.Elevation - vcen) <= GRID.ModelDZ_m_/2.));
            korners = [5,8,7,6];
    end
    numel(iplot)
    
    %% loop for animations
    %for k = 1:numel(ktimes)
    for k = 1
        nf=nf+1;figure;hold on;colormap(jet)
        
        if dodiffs == 1
            %% change in temperature T over a time interval           
            if upper(p_or_t) == 'P'
                ptlab = '\Delta P [Pa]';
            else
                ptlab = '\Delta T [degC]';
            end
            
            title(sprintf('%s over %d days from %s to %s \nX = %.1f Y = %.1f Z = %1.f E = %.1f N = %.1f V = %.1f [km]\n' ...
                ,ptlab ...
                ,days(PTFIELD.t(ktimes(2))-PTFIELD.t(ktimes(1)))...
                ,char(PTFIELD.t(ktimes(1))), char(PTFIELD.t(ktimes(2)))...
                ,xcen/1.e3,ycen/1.e3,zcen/1.e3 ...
                ,ecen/1.e3,ncen/1.e3,vcen/1.e3));
            if upper(p_or_t) == 'P'
                %this loop is slow
                for i=1:numel(iplot)
                    patch(GRID.CornersX(iplot(i),korners)/1000.,GRID.CornersY(iplot(i),korners)/1000.,GRID.CornersElevation(iplot(i),korners)/1000. ...
                        ,(PTFIELD.P(iplot(i),ktimes(2)) - PTFIELD.P(iplot(i),ktimes(1))) *ones(numel(korners),1));
                end
            else
                %this loop is slow
                for i=1:numel(iplot)
                    patch(GRID.CornersX(iplot(i),korners)/1000.,GRID.CornersY(iplot(i),korners)/1000.,GRID.CornersElevation(iplot(i),korners)/1000. ...
                        ,(PTFIELD.T(iplot(i),ktimes(2)) - PTFIELD.T(iplot(i),ktimes(1))) *ones(numel(korners),1));
                end
            end
        else
            %%  P or T at a single epoch
            if upper(p_or_t) == 'P'
                ptlab = 'P [Pa]';
            else
                ptlab = 'T [degC]';
            end
            
            title(sprintf('%s at %s X = %.1f Y = %.1f Z = %1.f E = %.1f N = %.1f V = %.1f [km]\n' ...
                ,ptlab ...
                ,char(PTFIELD.t(ktimes(k)))...
                ,xcen/1.e3,ycen/1.e3,zcen/1.e3 ...
                ,ecen/1.e3,ncen/1.e3,vcen/1.e3));
                       
            if upper(p_or_t) == 'P'
                %this loop is slow
                for i=1:numel(iplot)
                    patch(GRID.CornersX(iplot(i),korners)/1000.,GRID.CornersY(iplot(i),korners)/1000.,GRID.CornersElevation(iplot(i),korners)/1000. ...
                        ,PTFIELD.P(iplot(i),ktimes(k))*ones(numel(korners),1));
                end
            else
                %this loop is slow
                for i=1:numel(iplot)
                    patch(GRID.CornersX(iplot(i),korners)/1000.,GRID.CornersY(iplot(i),korners)/1000.,GRID.CornersElevation(iplot(i),korners)/1000. ...
                        ,PTFIELD.T(iplot(i),ktimes(k))*ones(numel(korners),1));
                end
            end
            %% try to unroll the loop
            %patch(GRID.CornersX(iplot,korners)/1000.,GRID.CornersY(iplot,korners)/1000.,GRID.CornersElevation(iplot,korners)/1000.,repmat(PTFIELD.T(iplot,ktimes(k)),1,numel(korners)));
            
        end
    end
    
    xlabel('ModelX [km]');
    ylabel('ModelY [km]');
    zlabel('Elevation [km]');
    
    switch j
        case 1
            view([-1,0,0]);
        case 2
            view([0,-1,0]);
        case 3
            view([0,0,1])
    end
    axis tight; axis equal;
    colorbar;
    
    %printjpg(sprintf('%s_View%1d_K%03d.jpg',mfilename,j,k));
    printpdf(sprintf('%s_View%1d_K%03d.pdf',mfilename,j,k));
end

return



