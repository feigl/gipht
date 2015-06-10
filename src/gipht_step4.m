% GIPHT_step4: make plots and images of gradients for sub-region
% 2014-06-03
fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

clearvars;
load


% only makes sense to do this for gradient from list
if ~ismember(pselect,[7,9])
    fprintf(1,'Skipping %s because pselect = %d\n',mfilename,pselect);
    return
end

% coordinates in meters
xmin=min(DST.x);
xmax=max(DST.x);
ymin=min(DST.y);
ymax=max(DST.y);
wesn = [xmin xmax ymin ymax];

phao = DST.phaobs;
% pass nulls to observations
if bitget(figopt,1) == 1
    % DEM elevation is zero where we have no data
    inull=find(abs(DST.z(i1:i2)) < 1);
    phao(inull) = NaN;
end

% residual gradients, without wrapping
gre0 = phao - mdl0;
gre1 = phao - mdl1;

% absolute deviations, without wrapping
gde0 = abs(gre0);
gde1 = abs(gre1);


% initialize symbols
for ii=1:8
    mysyms{ii} = '';
    marksizes(ii) = 1;
end
for ii=[5 6 7]
    mysyms{ii} = 'ko-';
    marksizes(ii) = 2;
end

% need to do this once
nlevels=64;cmap='jet';
tstr='East component of phase gradient after quad-tree partitioning (cycles/pixel)';
method=3;
zeroisnan=1;
%[imout, cmap, cuts]=histeq2(double(grx/DNPC),nlevels,method,cmap,tstr,zeroisnan);
imout = double(grx/DNPC);
cmap=colormap(cmap);
%cuts = nanmin(nanmin(imout)):nanmax(nanmax(imout))

% determine cut points from observed values - linear spacing
% nlevels = 64;
% cuts=linspace(nanmin(double(colvec(grx)/DNPC)),nanmax(double(colvec(grx)/DNPC)),nlevels);


% loop over pairs
for i = 1:np
    fprintf(1, 'Making gradient figures for pair %03d\n',i);
    
    % select row corresponding to this pair
    DD1 = DD(i,:);
    
    jmast =  abs(imast(i));
    jslav =  abs(islav(i));
    
    kmast = find(DD(i,:) == -1);
    kslav = find(DD(i,:) == +1);
    
    % indices into samples
    % New scheme using structures 2010-11-08
    kpairs = find(DST.k == i);
    i1 = DST.i(kpairs(1));
    i2 = DST.i(kpairs(end));
    %fprintf(1,'Starting stopping indices i1,i2 %d %d\n',i1,i2);
    
    titlestr = sprintf('Pair %3d orbs %5d %5d years %6.1f to %6.1f Dt = %.4f yr \nCost0 = %6.4f Cost1=%6.4f %s\n'...
        ,i,iuniqorbs(kmast),iuniqorbs(kslav)...
        ,tepochs(kmast),tepochs(kslav),tepochs(kslav)-tepochs(kmast)...
        ,nanmean(colvec(costs0(i1:i2)))/DNPC,nanmean(colvec(costs1(i1:i2)))/DNPC...
        ,strrep(runname,'_','\_'));
    
    % Make profiles %%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(1,'Making E-W profile for Pair %03d\n',i);
    nf=nf+1;h(nf)=figure;
    subplot(3,1,1);
    axis off;
    title(strcat(sprintf('Northing = %.3f km :',ycenter/1000),titlestr));

    subplot(3,1,2);
    % set up temporary plotting variables for observed values
    xt = DST.x(i1:i2)/1000; % easting in km
    %yt = phao(i1:i2)/DNPC;  % observed gradient in cycles per pixel
    yt = (mpercy/DNPC) * phao(i1:i2) ./ DST.dx(i1:i2);  % observed gradient in dimensionless units
    iok = find(isfinite(yt)==1);
    iok = intersect(iok,find(abs(yt)>0.0));
    iok = intersect(iok,find(abs(DST.y(i1:i2)-ycenter) < abs(5.0*DST.dy(i1:i2)))); % take pixels with 5 pixels widths of profile line
    
    if numel(iok) > 10
        %[dumb,iok] = sort(xt(iok));
        xt = xt(iok);
        yt = yt(iok);
        plot(xt,yt,'ro','MarkerFaceColor','r','MarkerSize',4);
        axis ij; % make range change increase downward
        hold on;
        % set up temporary plotting variables for modeled values
        xt = DST.x(i1:i2)/1000; % easting in km
        %yt = wrm1(i1:i2)/DNPC; % modeled gradient in cycles per pixel
        yt = (mpercy/DNPC) * wrm1(i1:i2) ./ DST.dx(i1:i2);  % modeled gradient in dimensionless units
        xt = xt(iok);
        yt = yt(iok);
        % sort in order of increasing x to connect lines
        %[dumb,isort] = sort(xt);
        % draw an envelope
        try
            isort = convhull(xt,yt);
        catch
            warning('convhull failed');
            isort = 1:numel(xt);
        end
        xt = xt(isort);
        yt = yt(isort);
        plot(xt,yt,'k-','LineWidth',1);
        set(gca,'FontName','Helvetica-Bold','Fontsize',12);
        legend('Observed','Modeled','Location','NorthEast');
        fixlabels('Easting (km)','%4.0f','gradient','');
        %title(strcat(sprintf('Northing = %.3f km :',ycenter/1000),titlestr));
        
        subplot(3,1,3);
        axis xy;
        % set up temporary plotting variables for residual values
        xt = DST.x(i1:i2)/1000; % easting in km
        %yt = wrm1(i1:i2)/DNPC;  % residual gradient in cycles per pixel
        yt = (mpercy/DNPC) * wrm1(i1:i2) ./ DST.dx(i1:i2);  % observed gradient in dimensionless units
        %iok = find(isfinite(yt)==1);
        %iok = intersect(iok,find(abs(DST.y(i1:i2)-ycenter) < abs(5*DST.dy(i1:i2))));
        %[dumb,iok] = sort(xt(iok));
        xt = xt(iok);
        yt = yt(iok);
        plot(xt,yt,'bo','MarkerFaceColor','b','MarkerSize',4);
        hold on;
        set(gca,'FontName','Helvetica-Bold','Fontsize',12);
        legend('Residual','Location','NorthEast');
        fixlabels('Easting (km)','%4.0f','gradient','');
        %title(strcat(sprintf('Northing = %.3f km :',ycenter/1000),titlestr));
        feval(printfun,sprintf('%s_PROFILEew_P%03d',runname,i));
        
    else
        warning(sprintf('Only %d points here. Not plotting.\n',numel(iok)));
    end
    
    
    % Make profiles %%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(1,'Making S-N profile for Pair %03d\n',i);
    nf=nf+1;h(nf)=figure;
    subplot(3,1,1);
    axis off;
    title(strcat(sprintf('Northing = %.3f km :',ycenter/1000),titlestr));

    subplot(3,1,2)
    % set up temporary plotting variables for observed values
    xt = DST.y(i1:i2)/1000; % northing in km
    %yt = phao(i1:i2)/DNPC;  % observed gradient in cycles per pixel
    yt = (mpercy/DNPC) * phao(i1:i2) ./ DST.dx(i1:i2);  % observed gradient in dimensionless units
    iok = find(isfinite(yt)==1);
    iok = intersect(iok,find(abs(yt)>0.0));
    iok = intersect(iok,find(abs(DST.x(i1:i2)-xcenter) < abs(5.0*DST.dx(i1:i2)))); % take pixels with 5 pixels widths of profile line
    if numel(iok) > 10
        
        %[dumb,iok] = sort(xt(iok));
        xt = xt(iok);
        yt = yt(iok);
        plot(xt,yt,'ro','MarkerFaceColor','r','MarkerSize',4);
        axis ij; % make range change increase downward
        hold on;
        % set up temporary plotting variables for modeled values
        xt = DST.y(i1:i2)/1000; % easting in km
        %yt = wrm1(i1:i2)/DNPC; % modeled gradient in cycles per pixel
        yt = (mpercy/DNPC) * wrm1(i1:i2) ./ DST.dx(i1:i2);  % modeled gradient in dimensionless units
        xt = xt(iok);
        yt = yt(iok);
        % sort in order of increasing x to connect lines
        %[dumb,isort] = sort(xt);
        % draw an envelope
        try
            isort = convhull(xt,yt);
        catch
            warning('convhull failed');
            isort = 1:numel(xt);
        end

        xt = xt(isort);
        yt = yt(isort);
        plot(xt,yt,'k-','LineWidth',1);
        set(gca,'FontName','Helvetica-Bold','Fontsize',12);
        legend('Observed','Modeled','Location','NorthEast');
        %title(strcat(sprintf('Easting = %.3f km :',xcenter/1000),titlestr));
        fixlabels('Northing (km)','%4.0f','Range gradient','');
        
        subplot(3,1,3);
        axis xy;
        % set up temporary plotting variables for residual values
        xt = DST.y(i1:i2)/1000; % easting in km
        %yt = wrm1(i1:i2)/DNPC;  % residual gradient in meters per pixel
        yt = mdl1(i1:i2)/DNPC;  % residual gradient in meters per pixel
        yt = (mpercy/DNPC) * wrm1(i1:i2) ./ DST.dx(i1:i2);  % residual gradient in dimensionless units
        %iok = find(isfinite(yt)==1);
        %iok = intersect(iok,find(abs(DST.y(i1:i2)-ycenter) < abs(5*DST.dy(i1:i2))));
        %[dumb,iok] = sort(xt(iok));
        xt = xt(iok);
        yt = yt(iok);
        plot(xt,yt,'bo','MarkerFaceColor','b','MarkerSize',4);
        hold on;
        set(gca,'FontName','Helvetica-Bold','Fontsize',12);
        legend('Residual','Location','NorthEast');
        fixlabels('Northing (km)','%4.0f','Range gradient','');
        %title(strcat(sprintf('Northing = %.3f km :',ycenter/1000),titlestr));
        feval(printfun,sprintf('%s_PROFILEsn_P%03d',runname,i));
        
    else
        warning(sprintf('Only %d points here. Not plotting.\n',numel(iok)));
    end


    
    % %     %  determine cut points from observed values - exp spacing
    cuts = colvec(quantile(phao(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),linspace(0.0,1.0,nlevels)));
        
    % use cuts to create indexed images
    zstr='gradient [dimensionless]'; drawcolorbar=0;
    
    % observed values
    tlA = 'Initial';
    tlE = 'Final';
    %imA=plot_patch_quads(phao(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imA=plot_patch_quads(phao(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imE=imA;
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlA,cuts(ceil(min(colvec(imA)))),cuts(ceil(max(colvec(imA)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlE,cuts(ceil(min(colvec(imE)))),cuts(ceil(max(colvec(imE)))));
    
%     % no need to wrap gradients 20120504
    tlB = 'Mdl0';
    tlF = 'Mdl1';
%     imB=plot_patch_quads(mdl0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imF=plot_patch_quads(mdl1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imB=plot_patch_quads(mdl0(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imF=plot_patch_quads(mdl1(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlB,cuts(ceil(min(colvec(imB)))),cuts(ceil(max(colvec(imB)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlF,cuts(ceil(min(colvec(imF)))),cuts(ceil(max(colvec(imF)))));
%    modeled values - wrapped
%     tlB = 'Mod0';
%     tlF = 'Mod1';
%     imB=plot_patch_quads(wrm0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imF=plot_patch_quads(wrm1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlB,cuts(ceil(min(colvec(imB)))),cuts(ceil(max(colvec(imB)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlF,cuts(ceil(min(colvec(imF)))),cuts(ceil(max(colvec(imF)))));
    
    %     % residual values, not wrapped
    tlC = 'Gre0';
    tlG = 'Gre1';
%     imC=plot_patch_quads(gre0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imG=plot_patch_quads(gre1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imC=plot_patch_quads(gre0(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imG=plot_patch_quads(gre1(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlC,cuts(ceil(min(colvec(imC)))),cuts(ceil(max(colvec(imC)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlG,cuts(ceil(min(colvec(imG)))),cuts(ceil(max(colvec(imG)))));
%     % residual values, with wrapping
%     tlC = 'Res0';
%     tlG = 'Res1';
%     imC=plot_patch_quads(res0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imG=plot_patch_quads(res1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlC,cuts(ceil(min(colvec(imC)))),cuts(ceil(max(colvec(imC)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlG,cuts(ceil(min(colvec(imG)))),cuts(ceil(max(colvec(imG)))));
    
%     % deviations, with abs function
%     tlD = 'Gde0';
%     tlH = 'Gde1';
%     imD=plot_patch_quads(gde0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imH=plot_patch_quads(gde1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlD,cuts(ceil(min(colvec(imD)))),cuts(ceil(max(colvec(imD)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlH,cuts(ceil(min(colvec(imH)))),cuts(ceil(max(colvec(imH)))));
    % cost values, with arc function
    tlD = 'Dev0';
    tlH = 'Dev1';
%     imD=plot_patch_quads(costs0(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     imH=plot_patch_quads(costs1(i1:i2)/DNPC,qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imD=plot_patch_quads(costs0(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
    imH=plot_patch_quads(costs1(i1:i2)*(mpercy/DNPC)./DST.dx(i1:i2),qii1(i1:i2),qii2(i1:i2),qjj1(i1:i2),qjj2(i1:i2),cuts,cmap,titlestr,zstr,drawcolorbar,nrsub,ncsub);
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlD,cuts(ceil(min(colvec(imD)))),cuts(ceil(max(colvec(imD)))));
%     fprintf(1,'%10s Min %12.6f Max %12.6f\n',tlH,cuts(ceil(min(colvec(imH)))),cuts(ceil(max(colvec(imH)))));
    
%     %Pass nulls from observations to residuals and totcosts
%     if bitget(figopt,1) == 1
% %         Pass nulls from observations
%         inull = [];
% %        inull = find(imA==1); % images contain indices to color map
% %         inull=union(inull,find(abs(imA) < 1/256.0));
% %         inull=union(inull,find(isfinite(imA)==0));
% %         inull=union(find(abs(imE) < 1/256.0),find(isfinite(imE)==0)); %
% %         areas with elevation less than 1 m
%         inull=union(inull,find(abs(DST.z(i1:i2)) < 1));
%         imA(inull) = NaN;
%         imB(inull) = NaN;
%         imC(inull) = NaN;
%         imD(inull) = NaN;
%         imE(inull) = NaN;
%         imF(inull) = NaN;
%         imG(inull) = NaN;
%         imH(inull) = NaN;
%     end

%     [imask2,jmask2] = find(amp > 10);
%     K = convhull(imask2,jmask2);    
%     [jmesh,imesh] = meshgrid(1:ncols,1:nrows);
%     imaskq = find(inpolygon(colvec(jmesh),colvec(imesh),jmask2(K),imask2(K)));
%     qsp(imaskq) = psp(imaskq);

    
    % 20131106
    if dl > 0
        imA = flipud(imA);
        imB = flipud(imB);
        imC = flipud(imC);
        imD = flipud(imD);
        imE = flipud(imE);
        imF = flipud(imF);
        imG = flipud(imG);
        imH = flipud(imH);
    end

% 20131120 Not needed for gradients  
%     % store images for big panel plot
%     oims(i,:,:) = imA;
%     mims(i,:,:) = imF;
%     rims(i,:,:) = imG;
%     cims(i,:,:) = imH;
%     tls{i} = titlestr;
    
    if pselect == 7 % flexible scale for gradient
        climit=cuts;
        ctab = cmapblackzero(1); % black at zero at bottom of color bar
        
        %hunt black patch bug 201204504
        %climit = [-1,0,+1];
        %ctab = colormap('jet');
    else
        climit=[-0.5, +0.5];
        ctab = cmapblackzero; % black at zero in central value
    end
    
    
    % make 8-panel plot in portrait
    nf=nf+1; h(nf)=utmimage8(imA,imB,imC,imD,imE,imF,imG,imH...
        ,tlA,tlB,tlC,tlD,tlE,tlF,tlG,tlH...
        ,wesn,strrep(titlestr,'\n',' '),climit,NaN,NaN,ctab,1,mysyms,marksizes);
    feval(printfun,sprintf('%s_%03d_8QAN',runname,i));
    
    % make 8-panel plot in portrait
    nf=nf+1; h(nf)=utmimage8landscape(imA,imB,imC,imD,imE,imF,imG,imH...
        ,tlA,tlB,tlC,tlD,tlE,tlF,tlG,tlH...
        ,wesn,strrep(titlestr,'\n',' '),climit,NaN,NaN,ctab,1,mysyms,marksizes);
    feval(printfun,sprintf('%s_%03d_8QANLS',runname,i));
    
    
    % calculate some derived parameters
    iq=iq+1;
    q0(iq) = max(max(abs(DST.mpercy(i1:i2).*phao(i1:i2)/DNPC./DST.dx(i1:i2))));
    q1(iq) = max(max(abs(DST.mpercy(i1:i2).*phao(i1:i2)/DNPC./DST.dx(i1:i2))));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Pair_%05d_MaxAbsGradientOBS',i);
    iq=iq+1;
    q0(iq) = max(max(abs(DST.mpercy(i1:i2).*mdl0(i1:i2)/DNPC./DST.dx(i1:i2))));
    q1(iq) = max(max(abs(DST.mpercy(i1:i2).*mdl1(i1:i2)/DNPC./DST.dx(i1:i2))));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Pair_%05d_MaxAbsGradientMOD',i);
    iq=iq+1;
    q0(iq) = max(max(abs(DST.mpercy(i1:i2).*res0(i1:i2)/DNPC./DST.dx(i1:i2))));
    q1(iq) = max(max(abs(DST.mpercy(i1:i2).*res1(i1:i2)/DNPC./DST.dx(i1:i2))));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Pair_%05d_MaxAbsGradientRES',i);
    iq=iq+1;
    q0(iq) = max(max(abs(DST.mpercy(i1:i2).*costs0(i1:i2)/DNPC./DST.dx(i1:i2))));
    q1(iq) = max(max(abs(DST.mpercy(i1:i2).*costs1(i1:i2)/DNPC./DST.dx(i1:i2))));
    qsig(iq) = NaN;
    qnames{iq} = sprintf('Pair_%05d_MaxAbsGradientDEV',i);
    
    
    %  PRINT OUT THE DERIVED PARAMETERS
    qnames = truncate_parameter_names(qnames);
    iq1 = iq2+1;
    iq2 = iq;
    uqb(iq1:iq2)=NaN;
    lqb(iq1:iq2)=NaN;
    qsig(iq1:iq2)=NaN;
    for j=iq1:iq2
        qflags{j} = 'D#';
        adj = q1(j)-q0(j);
        sadj = NaN;
        
        outfmt = getfmt(q1(j),qnames{j});
        
        fprintf(1        ,outfmt,qflags{i},j+mparam,qnames{j} ,q0(j),q1(j),adj,qsig(j),sadj,(uqb(j)-lqb(j))/2.0);
        fprintf(fidtxtout,outfmt,qflags{i},j+mparam,qnames{j}, q0(j),q1(j),adj,qsig(j),sadj,(uqb(j)-lqb(j))/2.0);
    end
end

clear phao imA imB imC imD imE imF imG imH;
clear h;

save('step4.mat');
save('qsave.mat','iq','iq1','iq2','qflags','qnames','q0','q1','qsig');

fprintf(1,'\n\n----------------   %s ended normally at %s ----------\n',upper(mfilename),datestr(now,31));

return



