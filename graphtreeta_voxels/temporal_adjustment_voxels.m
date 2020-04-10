%% temporal adjustment for Brady by voxel
% last updated 20191013 Elena Reinisch

%% initialize
clear all; close all;

child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))
verbose = 1;
nf=0;

load('VoxelMap.mat')

current_dir = pwd;
pdir = sprintf('%s/%s', pwd, 'Plots'); %name of folder to save prints to, if want in current working directory, replace with pwd
addpath(genpath(current_dir)); %calls on subdirectories as well (i.e. 'functions', etc)

printcleanfig = 1; % 1 = for papers, 0 = for research
ylab = 'Volume Change'; yunits = 'm^3 \times 10^6'; scalefactor = 1.e6;

% start text file
MSE = zeros(1656, 1);
SSWR = zeros(1656, 1);

% define temporal function
tfuncstxt = 'nsegs'
% tfuncstxt = 'rate'


%% Run through voxels for TA
for voxel_ind = 1:1656
    voxel_ind
    
    datfile=sprintf('brady_voxel%d_mest.txt', voxel_ind)
    
    dataset = datfile;
    
    titlestring = sprintf('Thermal Volumetric Strain Rate');
    
    
    insar_bperp = readtable('brady_msf_bperp.txt', 'ReadVariableNames', 0);
    
    %% Read Data file
    [VALS,NAMES] = read_table(datfile,0,5,5);
    [ndat,ndummy] = size(VALS);  % number of pairs
    tm = cal2datetime(VALS(:,1));  % master epoch in decimal years
    ts = cal2datetime(VALS(:,2));  % slave epoch in decimal years
    
    ndat = numel(tm);
    sat_id = zeros(ndat, 1);
    [TSX_ind, mdum] = find(char(NAMES(:,3)) == (NAMES(end,3)));
    [ERS2_ind, mdum] = find(char(NAMES(:,3)) == char(NAMES(1, 3)));
    [ALOS_ind, mdum] = find(char(NAMES(:,3)) == char(NAMES(6, 3)));
    sat_id(TSX_ind) = 3;
    sat_id(ERS2_ind) = 1;
    sat_id(ALOS_ind) = 2;
    
    
    ndat = numel(tm);
    
    tm_dt = tm;
    tm_dt.Format = 'yyy-MM-dd';
    ts_dt = ts;
    ts_dt.Format = 'yyy-MM-dd';
    tm = dyear(year(tm), month(tm), day(tm));
    ts = dyear(year(ts), month(ts), day(ts));
    
    Qrate = VALS(:,3);
    rs = (abs(VALS(:,4)))*sqrt(1e11);% scale based on original MSE
       
    % unique epochs in years
    tu=sort(unique([tm ts]));
    tu_dt = sort(unique([tm_dt ts_dt]));
    tu_dt.Format = 'yyy-MM-dd';
    
    tspan = sprintf('%10.4f to %10.4f\n',min(tu),max(tu));
    
    % midpoints of each interval
    tmid = (tm+ts)/2;
    
    % half interval
    th = (ts-tm)/2.;
    
    % initial epoch
    t0 = min(tu);

    
    %% Plot rates
    if voxel_ind == 1
        h_fig = figure;
        set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
        set(h_fig,'DefaultTextInterpreter','tex');
        hold on;

        indother = 1:numel(tm);
        
        pest1 = mean(Qrate/scalefactor);
        psig1 = std(Qrate/scalefactor);
        
        errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
        errorbar_plus2(tmid(indother),Qrate(indother)/scalefactor,th(indother),rs(indother)/scalefactor,'gs',5);
       ymin = nanmin(Qrate-rs)/scalefactor;
        ymax = nanmax(Qrate+rs)/scalefactor;
       errorbar_plus2((max(tu)+min(tu))/2.,pest1,(max(tu)-min(tu))/2.,psig1,'ko-',5);
        
        % draw tbreaks in green
        %tbreaks = [tu(1), dyear(2016, 02, 01), dyear(2016, 03, 14), dyear(2016, 03, 24), dyear(2016, 08, 23), tu(end)];
        tbreaks = [tu(1), dyear(2017, 04, 27), tu(end)];
        
        for i=1:numel(tbreaks)
            plot([tbreaks(i) tbreaks(i)],[ymin ymax],'k--');
        end
        hold off;
        xlabel('Date');
        ylabel(sprintf('Rate of %s [%s/year]',ylab,yunits));
        set(gca, 'XTickLabel', char(tu_dt), 'XTickLabelRotation', -75)
        if printcleanfig == 0
            title(strcat(titlestring,sprintf('\n Mean rate is %#10.4e +/- %#10.4e [%s/year]',pest1(1)/scalefactor,psig1(1)/scalefactor,yunits)));
            
            printpdf(sprintf('%s_%sRATES.pdf',mfilename,dataset));
        elseif printcleanfig == 1
            set(gca, 'FontSize', 14)
            title('Estimated Rates of Volume Change for TSX pairs from 2017')
            [axy, axm, axd] = dyear2yyyymmdd(get(gca, 'XTick'));
            tick_labels = datetime(axy, axm, axd);
            set(gca, 'XTickLabel',datestr(tick_labels, 'yyyy-mmm-dd'), 'XTickLabelRotation', -75) %char(tick_labels)
            printeps_e('Coso_mogi_RATES.eps');
            
        end
    end
    
    %% Convert rates into differences over time intervals
    % this converts the volume from rate (m^3/year) to differential (m^3)
    Qdiff = Qrate.* (ts - tm); % Qdiff is differential quantity observed over time interval
    Qdsig = rs.* (ts - tm); % Qdsig is uncertainty (1-sigma) for differential observation
    
    %% find and plot trees
    if voxel_ind == 1
        dispflag = 0;
        
        [tm_y, tm_m, tm_d] = dyear2date(tm);
        [ts_y, ts_m, ts_d] = dyear2date(ts);
        tm_dstruct = datetime(tm_y, tm_m, tm_d);
        ts_dstruct = datetime(ts_y, ts_m, ts_d);
        
        [trees, DD tepochs, iepochs, iuniqorbs, uniqdates] = findtrees2(tm_dstruct',ts_dstruct');
        
        [ndummy,mepochs] = size(DD);
        fprintf(1,'Number of unique epochs %d\n',mepochs);
        [ntrees,ndummy] = size(trees);
        fprintf(1,'Number of distinct trees %d\n',ntrees);
        
        %% naive temporal adjustment
        Qdiff2 = adjustbp(tepochs,DD, insar_bperp.Var3, trees,iuniqorbs, uniqdates);
        cal_date = [month(tepochs) day(tepochs)];
        % plot
        
        tepochs = dyear(year(tepochs), month(tepochs), day(tepochs));
        
        if printcleanfig == 0 %ktours = plotbp...
            plotbp(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), cal_date);
            printpdf(sprintf('%s_%s_TREES.pdf',mfilename,dataset),h);
        elseif printcleanfig == 1
            figure;
            set(gca,'FontName','Helvetica','Fontweight','Bold','FontSize',12);
            set(gcf,'DefaultTextInterpreter','tex');
            plotbp_bradyMSF(tepochs, Qdiff2, DD, trees, iuniqorbs, uniqdates, 0,sprintf('%s [%s]',ylab,yunits), sat_id);
            ylabel('orbital separation [m]')
            h =  legend('ALOS', 'TSX', 'ERS2')
            h.FontSize = 16;

        end
    end
    
    %% loop over temporal fitting functions
    switch tfuncstxt
        case 'rate'
            tfuncs = {'rate'};
        case 'nsegs'
            tfuncs = {'nsegs'};
    end

    % dimension arrays
    mses = nan(numel(tfuncs),1);
    mparams = nan(numel(tfuncs),1);
    
    for k = 1:numel(tfuncs)
        tfunc = tfuncs{k};
        
        fprintf(1,'\n-----------\n');
        
        fprintf(1,'Time function is now %s\n',tfunc);
        
        
        switch tfunc
            case {'nsegs','step-pwl','step-sec', 'nsegs0'};
                metaparams = nan;
                tbreaks = [tu(1), dyear([2004:2018],1*ones(size([2004:2018])), 1*ones(size([2004:2018]))), tu(end)];
                tbreaks = unique(sort(tbreaks));
            otherwise
                metaparams = nan;
                % few breaks
                tbreaks = [];
                tbreaks(end+1) = min(tu);
                tbreaks(end+1) = max(tu);
        end
        
        %% perform temporal adjustment
        tbreaks = colvec(sort(unique(tbreaks)))
       [pest, psig, mse, Qdmod, tfit, pfit, sigl, sigu, rd, V, G, sswr, Vx, var, res_n] = temporal_adjustment(Qdiff,Qdsig,tm,ts,tbreaks,tfunc,metaparams);
        
        mparam = numel(pest);
        %%
        % store statistics
        mses(k) = mse;
        mparams(k) = mparam;

    end
    
    % predict estimates for each epoch
    Qdmod_tu = [];
    Qdsigl_tu = [];
    Qdsigu_tu = [];
    for i=1:numel(tu)
        tf = time_function(tfunc, tu(i),   tbreaks, metaparams); % find time value for given parameterization
        dm=0;
        dsl = 0;
        dsu = 0;
        for j=1:mparams
            dm  = dm   + tf(j) * pest(j); % interpolate to find displacement per epoch
            dsl = dsl + tf(j) * (pest(j)-psig(j));
            dsu = dsu + tf(j) * (pest(j)+psig(j));
        end
        Qdmod_tu(i,1) = dm; % estimated model displacement
        Qdsigl_tu(i,1) = dsl; % lower 1-sigma uncerainty bounds
        Qdsigu_tu(i,1) = dsu; % upper 1-sigma uncertainty bounds
    end % loop over epochs
    
    tu_cal = dyear2caldate(tu);
    tu_cal = year(tu_dt)*1e4 + month(tu_dt)*1e2 + day(tu_dt);
    
    % save mse
    MSE(voxel_ind) = mse;
    SSWR(voxel_ind) = sswr;
    
    % diary off
    dlmwrite(char(strcat('TA_cube_results_', tfuncstxt, '.txt')), [voxel_ind*ones(size(tu_cal)), tu_cal, Qdmod_tu, Qdsigl_tu, Qdsigu_tu, VoxelMap.centx_utm(voxel_ind)*ones(size(tu_cal)), VoxelMap.centy_utm(voxel_ind)*ones(size(tu_cal))], '-append', 'precision','%.12f')
    
end

dlmwrite(char(strcat('TA_cube_mse_', tfuncstxt, '.txt')), [[1:1656]', MSE], 'precision', '%5.5e')
dlmwrite(char(strcat('TA_cube_sswr_', tfuncstxt, '.txt')), [[1:1656]', SSWR], 'precision', '%5.5e')

return

