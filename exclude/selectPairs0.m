function pairListFilenameOutput = selectPairs(pairListFileNameInput,dataSet,doPlots,CRITERIA,DATES)
%function [] = avg_rate_range_grds( pair_list_file, metersperfringe )

% initialize
nf = 0; % number of figures
criteriaFieldNames = {'maxBperp','minDays','maxDays','dirName','fileNameEnding'}
datesFieldNames = {'yyyymmdd0','yyyymmdd1','yyyymmdd2'}

narginchk(1,5);

if nargin < 2
    dataSet = 'ALL';
end

if nargin < 3
    doPlots = 1;
end

if nargin < 4
    CRITERIA = struct();
end

if nargin < 5
    DATES = struct();
end

%% unpack Criteria
for i=1:numel(criteriaFieldNames)
    criteriaFieldName1 = criteriaFieldNames{i};
    if isfield(CRITERIA,criteriaFieldName1) == false
        switch criteriaFieldName1
            case 'maxBperp'
                CRITERIA.(criteriaFieldName1) = 100; % max orbital separation in meters
            case 'minDays'
                CRITERIA.(criteriaFieldName1) = 1; % mininum time span in days
            case 'maxDays'
                CRITERIA.(criteriaFieldName1) = 10*365; % maximum time span in days
            case 'dirName'
                CRITERIA.(criteriaFieldName1) = '../intf';
            case 'MSF'
                CRITERIA.(criteriaFieldName1) = true;
            case 'fileNameEnding'
                %fileNameEnding = 'phasefilt_mask_utm.grd'; % wrapped phase in radians GIPhT idatatype = 0
                %fileNameEnding = 'drhomaskd_utm.grd';      % range change in meters GIPhT idatatype = 1
                %fileNameEnding = 'dgrad_utm.grd';          % range gradient GIPhT idatatype = 1   does not work yet
                CRITERIA.(criteriaFieldName1) = 'drhomaskd_utm.grd';      % range change in meters GIPhT idatatype = 1
            otherwise
                error(sprintf('Unknown case %s for criteria\n',criteriaFieldName1));
        end
    end
end
CRITERIA




%% read file produced by egenerate_pairlist.sh
%
% Script for generating master pairlist for sat/track/site.
% Args: egenerate_pairlist.sh [site] [sat] [trk] [ALOS optional - frame number]
% Example:
% egenerate_pairlist.sh brady TSX T53
% Example: ALOS from frame 1190
% egenerate_pairlist.sh malas ALOS T244 1190
% ! rsync -rav maule.ssec.wisc.edu:/s21/insar/TSX/T91/preproc/TSX_T91_sanem_pairs.txt .
%pairListFileName = 'TSX_T91_sanem_pairs.txt';
% 20200408 Kurt Feigl include pairs with  20200306 and 20200328
% ! rsync -rav maule.ssec.wisc.edu:/s21/insar/TSX/T30/preproc/TSX_T30_forge_pairs.txt .
%pairListFileNameInput = 'TSX_T30_forge_pairs.txt';

TpairsIn = readFlatFileAsTable(pairListFileNameInput,1);
%TpairsIn = readtable(pairListFileNameInput)
[nrows, ncols0] = size(TpairsIn);

% skip this if only plots are requested
if doPlots < 2
    %% add 2 columns to table for pair name and file names
    pairnames = cell(nrows,1);
    filenames = cell(nrows,1);
    % %dirname = '../intf'; % unwrapping failed for some of these
    % dirName = '../intf1';
    % idatatype = 1;
    % switch idatatype
    %     case 0
    %         fileNameEnding = 'phasefilt_mask_utm.grd';
    %     case 1
    %         fileNameEnding = 'drhomaskd_utm.grd';
    %     case 2
    %         fileNameEnding = 'dgrad_utm.grd' % does not work yet
    %     otherwise
    %         error(sprintf('unknown idatatype = %d\n',idatatype));
    % end
    
    for i=1:nrows
        mast1 = sprintf('%8d',TpairsIn.mast(i));
        slav1 = sprintf('%8d',TpairsIn.slav(i));
        pairname1 = sprintf('%s_%8s',mast1,slav1);
        filename1 = strcat(CRITERIA.dirName,filesep,sprintf('In%s',pairname1),filesep,CRITERIA.fileNameEnding);
        pairnames{i} = pairname1;
        filenames{i} = filename1;
    end
    TpairsIn = [TpairsIn, table(pairnames,'VariableNames', {'pairname'})];
    TpairsIn = [TpairsIn, table(filenames,'VariableNames', {'filename'})];
    
    %% Add columns for dates
    % datetime structures
    TpairsIn = [TpairsIn, table(cal2datetime(TpairsIn.mast),'VariableNames', {'tm_dstruct'})];
    TpairsIn = [TpairsIn, table(cal2datetime(TpairsIn.slav),'VariableNames', {'ts_dstruct'})];
    
    % decimal years
    tm_dyear = dyear(year(TpairsIn.tm_dstruct),      month(TpairsIn.tm_dstruct),      day(TpairsIn.tm_dstruct));
    ts_dyear = dyear(year(TpairsIn.ts_dstruct),      month(TpairsIn.ts_dstruct),      day(TpairsIn.ts_dstruct));
    
    TpairsIn = [TpairsIn, table(tm_dyear,'VariableNames', {'tm_dyear'})];
    TpairsIn = [TpairsIn, table(ts_dyear,'VariableNames', {'ts_dyear'})];
    clear tm_dyear ts_dyear;   
end


%% find unique
[TpairsIn,ia,ic] = unique(TpairsIn,'rows');

% find duplicates
id = setdiff(ia,ic);
if numel(id) > 0
    fprintf(1,'WARNING: found and deleted duplicate rows:\n');
    TpairsIn(id,:);
end
TpairsIn

% count
npairs = numel(table2array(TpairsIn(:,1)))

%% select dates

% unpack Dates
for i=1:numel(datesFieldNames)
    datesFieldName1 = datesFieldNames{i}
    if isfield(DATES,datesFieldName1) == false
        switch i
            case 1  % origin epoch
                % round to nearest year
                yyyymmdd0 = nanmin([TpairsIn.mast, TpairsIn.slav]);
                yyyy=floor(yyyymmdd0/10000)*1000.
                yyyymmdd0 = yyyy + 0101;
                DATES.(datesFieldName1) = yyyymmdd0;
            case 2  % first acquisition date for stack
                yyyymmdd1 = nanmin([TpairsIn.mast, TpairsIn.slav]);
                DATES.(datesFieldName1) = yyyymmdd1;
            case 3  % last acquisition date for stack
                yyyymmdd2 = nanmax([TpairsIn.mast, TpairsIn.slav]);
                DATES.(datesFieldName1) = yyyymmdd2;
            otherwise
                error(sprintf('Unknown case i %d for date\n',i));
        end
    end
end


% timeSpans = {'ALL','FY20Q2'};
% timeSpans = timeSpans(1); % not enough data for early
%for itspan = 1:numel(timeSpans)
%close all;
timeSpan1 = sprintf('%8dTo%8d',DATES.yyyymmdd1,DATES.yyyymmdd2);

% switch timeSpan1
%     case 'ALL' %
%         yyyymmdd1 = 20161108; % min(Tpairs0.mast)
%         yyyymmdd2 = max(Tpairs0.mast) %  20200324
%     case 'FY20Q2' % FY'20 Q2 (ending March 30, 2020)
%         % To do so, we have analyzed the SAR data from early January 2019 (20190131) through March 2020 (20200324).
%         yyyymmdd1 = 20190131;
%         yyyymmdd2 = 20200324; %
%     otherwise
%         error(sprintf('unknown timeSpan1 %s\n',timeSpan1));
% end

fprintf(1,'Pruning dates between %d and %d\n',DATES.yyyymmdd1,DATES.yyyymmdd2);


%% work with logical variables
isok1 = true(npairs,1);
isok1 = and(isok1,isfinite(TpairsIn.mast));
isok1 = and(isok1,isfinite(TpairsIn.slav));
isok1 = and(isok1,(TpairsIn.mast >  DATES.yyyymmdd1));   % master date
isok1 = and(isok1,(TpairsIn.slav >  DATES.yyyymmdd1));   % slave  date
isok1 = and(isok1,(TpairsIn.mast <= DATES.yyyymmdd2));   % master date
isok1 = and(isok1,(TpairsIn.slav <= DATES.yyyymmdd2));   % slave  date

% convert to indices
iok1 = find(isok1 == true);


%% check that file exists
iok2 = zeros(npairs,1);
fnames = TpairsIn.filename;
kount = 0;
for i=1:npairs
    fname1 = char(fnames{i});
    if exist(fname1,'file') == 2
        kount=kount+1;
        iok2(kount) = i;
    else
        fprintf(1,'WARNING: could not open grd file named fname1 = %s . Skipping.\n',fname1);
    end
end

%% prune
iok = intersect(iok1,iok2);
TpairsOut = TpairsIn(iok,:);
[npairs, ncols] = size(TpairsOut);

%% Analyze graphs for several subsets of the data
%dataSets = {'ALL','MSF'};
%dataSets = {'ALL'};
%for idataSet = 1:numel(dataSets)
%
%dataSet1 = dataSets{idataSet}

if (npairs > 0) == false
    error(sprintf('Too few pairs npairs = %d\n',npairs));
else
    %% get unique list of dates, sorted
    tu_cal   = unique([TpairsOut.mast;  TpairsOut.slav]);  % unique epochs
    tu_dyear = unique([TpairsOut.tm_dyear, TpairsOut.ts_dyear]);
    tu_dstruct = unique([TpairsOut.tm_dstruct, TpairsOut.ts_dstruct]);
    
    ntu = numel(tu_cal);
    tu_str = cell(ntu,1);
    for i=1:ntu
        tu_str{i} = sprintf('%8d',tu_cal(i));
    end
    
    
    %% list of edges by their calendar date YYYYMMDD
    Edges = [TpairsOut.mast,TpairsOut.slav];
    [npairs0,ncols] = size(Edges);
    if ncols == 2
        npairs0
    else
        ncols
        error('ncols must equal 2.');
    end
    
    %% find edge-vertex incidence matrix
    [Qiev] = edges_to_incidence(Edges);
    [npairs1,nepochs] = size(Qiev);
    fprintf(1,'Number of unique epochs %d\n',nepochs);
    
    %% list of edges by their indices
    iedges=incidence_to_edges(Qiev);
    
    %% find list of trees indexing by vertices
    verbose = 0;
    itrees = find_trees_from_incidence(Qiev,verbose);
    
    [ntrees,ndummy] = size(itrees);
    fprintf(1,'Number of distinct trees %d\n',ntrees);
    
    %% find indices of pairs
    ipairs = trees_to_pairs(itrees,iedges);
    [npairs2,ncols] = size(ipairs);
    if ncols ~= 1
        ncols
        error('ncols must equal 1.');
    end
    
    %% sanity check
    if npairs1 == npairs2
        npairs = npairs1;
    else
        npairs1
        npairs2
        error('miscount');
    end
    
    
    
    %% choose criterion for adjustment
    %criteria = {'Bperp','BpDt'};
    %     criteria = {'BpDt'};
    %     for icriterion = 1:numel(criteria)
    criterion1 = 'BpDt'
    switch criterion1
        case 'Bperp'
            %Use absolute value of Bperp for temporal adjustment
            Weights = abs(TpairsOut.bperp);
            ylab = 'Bperp';
            yunits = 'm';
        case 'Dt'  % useless because graph is a straight line
            % Use time span in days
            Weights = abs(TpairsOut.dt);
            ylab = '\Delta t';
            yunits = 'days';
        case 'BpDt'
            %Use product of time span and Bperp
            %Weights = abs(zscore(Tpairs.bperp) .* zscore(Tpairs.dt));
            Weights = zscore(TpairsOut.bperp) .* zscore(TpairsOut.dt);
            ylab = 'normalized Bp * \Delta t';
            yunits = 'dimensionless';
            %                     case 'BpDtZs'
            %                         %Use product of time span and Bperp and Zscore of "stable" area
            %                         Weights = abs(zscore(Tpairs.bperp) .* zscore(Tpairs.dt) .* zscore(Tpairs.meanRate));
            %                         ylab = 'normalized Bp * \Delta t';
            %                         yunits = 'dimensionless';
        otherwise
            error(sprintf('unknown criterion %s\n', criterion1));
    end
    
    if doPlots
        
        %% make title for plots
        titlestr = sprintf('%s %s %s',strtrim(timeSpan1),strtrim(dataSet),strtrim(criterion1));
        
        %% estimate the Yvalues for the graph from the weights
        % this centers the yvalues
        YvalEst = invert_incidence_matrix(Weights,Qiev,itrees);
        
        %% make a graph using Matlab functions
        DGRAPH = digraph(iedges(:,1),iedges(:,2),Weights,tu_str);
        Nodes = table2array(DGRAPH.Nodes);
        
        %% calculate importance of nodes
        importance = centrality(DGRAPH,'authorities','Importance',abs(DGRAPH.Edges.Weight));
        fprintf(1,'i,Nodes{i},tu_cal(i),importance(i) \n');
        for i=1:nepochs
            fprintf(1,'%3d %8s %8d %10.4f\n',i,Nodes{i},tu_cal(i),importance(i));
        end
        
        %% is it a DAG
        idagness = isdag(DGRAPH);
        
        %% plot trees
        xlab = 'Date';
        xunits = 'year';
        xlabstr = sprintf('%s [%s]',xlab,xunits);
        ylabstr = sprintf('%s [%s]',ylab,yunits);
        
        
        % plot graph using Matlab function with labels
        %'auto', 'circle', 'force', 'layered', 'subspace', 'subspace3', 'force3'
        nf=nf+1;figure;
        %layout = 'force';
        layout = 'auto';
        %layout = 'circle';
        %layout = 'layered';
        %layout = 'subspace';
        plot(DGRAPH,'Layout',layout,'EdgeLabel',DGRAPH.Edges.Weight);
        title(titlestr);
        %                 printpdf(sprintf('%s_trees1_%s_%s_%s.pdf',mfilename ...
        %                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
        fname2 = 'digraph1';
        PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
        printpdf(PDFfileName);
        
        % plot graph using Matlab function with labels
        nf=nf+1;figure;
        plot(DGRAPH,'EdgeLabel',DGRAPH.Edges.Weight,'XData',tu_dyear, 'YData',YvalEst);
        xlabel(xlabstr);
        ylabel(ylabstr);
        title(titlestr);
        %                 printpdf(sprintf('%s_trees2_%s_%s_%s.pdf',mfilename ...
        %                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
        fname2 = 'digraph2';
        PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
        printpdf(PDFfileName);
        
        
        
        % plot graph using function from GraphTreeTA
        nf=nf+1;Hfig = plot_trees3(tu_dyear, YvalEst, Qiev, itrees, xlab, ylabstr, titlestr);
        %                 printpdf(sprintf('%s_trees3_%s_%s_%s.pdf',mfilename ...
        %                     ,strtrim(timeSpan1),strtrim(dataSet1),strtrim(criterion1)));
        fname2 = 'plot_trees3';
        PDFfileName = sprintf('%s_Fig%03d_%s.pdf',mfilename,nf,fname2);
        printpdf(PDFfileName);
        
    end
    
    %         %% write a table
    %         OutTableFileName = sprintf('%s_%s_%s_%s.ipairs',mfilename,strtrim(timeSpan1),strtrim(dataSet),strtrim(criterion1))
    %         %% TODO why not write MSF table?
    %         writetable(TpairsOut(1:npairs,1:ncols0),OutTableFileName,'FileType','text','Delimiter','space');
    
    %
    %% find mininum spanning forest
    %if contains(dataSet, 'MSF') == true
    if CRITERIA.MSF 
        rs = abs(Weights);
        scalefactor = 1;
        [i_msf, i_rep, nCycles, QievMSF, TreesMSF, WeightsMSF] = findMininumSpanningForest(Qiev, itrees, Weights, scalefactor);
        nCycles
        % prune to MSF
        TpairsOut = TpairsOut(i_msf,:);
    end
    
    
    [npairs, ncolsmax] = size(TpairsOut);
    
    %% check that arrays are OK
    % get information about first file
    %INFO=grdinfo3(char(grd_list(1)));
    % fname1 = TpairsOut.filename{1}
    % INFO1=grdinfo3(fname1);
    % [xgrd,ygrd,drho] = grdread3((fname1));
    % nf=nf+1;map_grd(fname1);
    
    
end

if doPlots == 2
    pairListFilenameOutput = NaN;
else
    pairListFilenameOutput = sprintf('%s_%s_%3s_ipairs.txt' ...
        ,strrep(pairListFileNameInput,'.csv','') ...
        ,dataSet ...
        ,'select');
    writetable(TpairsOut,pairListFilenameOutput);
    fprintf(1,'%s: successfully wrote new table named: %s\n',mfilename,pairListFilenameOutput);  
end


return
end






