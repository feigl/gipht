function [pdfFileNames] = scatter_map6(OBS,TMAP,do6panel,iQuantities,iosub,iPlot,tstr1,subset)
%% map quantities as scatter plots with pixels as circles in 6 panels
nf=0;


switch iPlot
    case 1  % maps
        pw = 0.25 % width for subplot in normalized units
        ph = 0.25 % height for subplot in normalized units
        bx = 0.03 % border between subplots in normalized units
        by = 0.15 % border between subplots in normalized units
    case 2  % histograms
        pw = 0.25 % width for subplot in normalized units
        ph = 0.25 % height for subplot in normalized units
        bx = 0.07 % border between subplots in normalized units
        by = 0.15 % border between subplots in normalized units
    otherwise
        error(sprintf('unknown iPlot %d\n',iPlot));
end


if do6panel == 1
    nf=nf+1;
    figure;
    %figure('position',[0 0 1 1]); % full screen
end

for iQuantity = iQuantities
    switch iQuantity
        case 1
            qstr = 'OBS';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = OBS.vlos;
            pos = [1*bx, 2*by+ph, pw, ph];
        case 2
            qstr = 'STD';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = OBS.qlos;
            pos = [1*bx, by, pw, ph];
        case 3
            qstr = 'MOD1';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = MOD1;
            pos = [2*bx+pw, 2*by+ph, pw, ph];
         case 4
            qstr = 'MOD2';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = MOD2;
            pos = [2*bx+pw, by, pw, ph];
        case 5
            qstr = 'RES1';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = RES1;
            pos = [3*bx+2*pw, 2*by+ph, pw, ph];
        case 6
            qstr = 'RES2';ustr='mm/year';vscl=1.e3;ylab='LOS rate';
            V = RES2;
            pos = [3*bx+2*pw, by, pw, ph];
        case 7 % same as MOD1
            qstr = 'TOPO';ustr='m';vscl=1.0;ylab='elevation';
            V = OBS.helv;
            pos = [2*bx+pw, 2*by+ph, pw, ph];
        otherwise
            error('unknown iQuantity');
    end
    
    if do6panel == 1
        %subplot(2,3,iplot);
        subplot('position',pos);
        tstr = sprintf('%s',qstr);
    else
        nf=nf+1;figure;
        tstr = sprintf('%s\n%s', qstr,tstr1);
    end
    
    switch iPlot
        case 1 % make map
            pstr = 'MAP1';
            circleSize = 10;
            scatter(OBS.eutm(iosub)/1.e3,OBS.nutm(iosub)/1.e3,circleSize,vscl*V(iosub),'filled'); % plot in mm/year
            hold on;
%             if iQuantity == 1
%                 plot(OBS.eutm(isamp)/1.e3,OBS.nutm(isamp)/1.e3,'k.'); % coordinates of pixels
%             end
%             if iQuantity == 3 || iQuantity == 4
%                 plot(VOXX.Xutm(ivsub)/1.e3,VOXX.Yutm(ivsub)/1.e3,'k+'); % coordinates of voxels
%                 plot(NODE.Xutm(insub)/1.e3,NODE.Yutm(insub)/1.e3,'ws'); % coordinates of nodes
%             end
%             plot(ResX/1.e3,ResY/1.e3,'wo','MarkerSize',5,'LineWidth',2); % center of modeled reservoir
            
            axis equal
            axis tight
            
            % plot wells
            if iQuantity == 1
                wellTypes = {'Production','Idle','Injection'};
                for i=1:numel(wellTypes)
                    wellType1=wellTypes{i}
                    switch wellType1
                        case 'Production'
                            sym1 = '^'; color1='r';
                        case 'Injection'
                            sym1 = 'v'; color1='b';
                        case 'Idle'
                            sym1 = 'o'; color1='g';
                        otherwise
                            warning(sprintf('Unknown wellType1 %s\n',wellType1));
                    end
                    ii=find(contains(TMAP.Utility,wellType1));
                    if numel(ii) > 0
                        plot(TMAP.Easting(ii)/1.e3,TMAP.Northing(ii)/1.e3,sprintf('k%s',sym1),'MarkerFaceColor',color1,'MarkerSize',10);
                    end
                end
            end
            
            
            if do6panel == 1 && iQuantity ~= 2
                axis off
            end
            
            % arrange colors
            colormap('jet');
            % set limits for color table
            if do6panel == 1
                % same for all panels
                maxabs=max([abs(nanmin(OBS.vlos(iosub))), abs(nanmax(OBS.vlos(iosub)))]);
            else
                maxabs=max([abs(nanmin(V(iosub))), abs(nanmax(V(iosub)))]);
            end

            if iQuantity < 7
                climits = 1.e3*[-1*maxabs,+1*maxabs]; % in mm
                caxis(climits); % clip color bar
            end
            
            % axis and title
            xlabel('UTM easting [km]');
            ylabel('UTM northing [km]');
            title(tstr,'Interpreter','none');
            
            if (do6panel == 1 && iQuantity == 6) || do6panel == 0
                CB=colorbar;
                %CB.Label.String = 'LOS rate';
                %title(CB,'mm/year');
                title(CB,ustr);
                %ylabel(CB,'LOS rate');
                ylabel(CB,ylab);
            end
            
            
        case 2 % make histogram
            if ismember(iQuantity,[1:6])
                % make histograms of residuals
                pstr = 'HIST';
                histogram(1.e3*V(isamp));
                axis tight;
                title(tstr,'Interpreter','none');
                xlabel('LOS velocity [mm/year]');
                ylabel('Count')
            end
        otherwise
            error(sprintf('unknown iPlot %d\n',iPlot));
    end
    
    if do6panel == 0
        pdfFileNames{nf} = sprintf('%s_%s_%s_%s.pdf',strrep(tstr1,' ','_'),pstr,qstr,subset);
        printpdf(pdfFileNames{nf});
    end
    
end
if do6panel == 1
    pdfFileNames{nf} = sprintf('%s_%s_%s_%s.pdf',strrep(tstr1,' ','_'),pstr,'6PAN',subset);
    %printpdf(sprintf('%s_%s_%s.pdf',runname,'6PAN',subset));
end
pdfFileNames=transpose(pdfFileNames); % column vector
return
end
