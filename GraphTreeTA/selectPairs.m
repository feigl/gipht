function Tpairs = selectPairs(Tpairs,dataSet,doPlots,CRITERIA,DATES)
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
for iPair=1:numel(criteriaFieldNames)
    criteriaFieldName1 = criteriaFieldNames{iPair};
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

% TpairsIn = readFlatFileAsTable(pairListFileNameInput,1);
% %TpairsIn = readtable(pairListFileNameInput)
[nrows, ncols0] = size(Tpairs);

% skip this if only plots are requested

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

for iPair=1:nrows
    mast1 = sprintf('%8d',Tpairs.mast(iPair));
    slav1 = sprintf('%8d',Tpairs.slav(iPair));
    pairname1 = sprintf('%s_%8s',mast1,slav1);
    filename1 = strcat(CRITERIA.dirName,filesep,sprintf('In%s',pairname1),filesep,CRITERIA.fileNameEnding);
    pairnames{iPair} = pairname1;
    filenames{iPair} = filename1;
end
Tpairs = [Tpairs, table(pairnames,'VariableNames', {'pairname'})];
Tpairs = [Tpairs, table(filenames,'VariableNames', {'filename'})];

%% Add columns for dates
% datetime structures
Tpairs = [Tpairs, table(cal2datetime(Tpairs.mast),'VariableNames', {'tm_dstruct'})];
Tpairs = [Tpairs, table(cal2datetime(Tpairs.slav),'VariableNames', {'ts_dstruct'})];

% decimal years
tm_dyear = dyear(year(Tpairs.tm_dstruct),      month(Tpairs.tm_dstruct),      day(Tpairs.tm_dstruct));
ts_dyear = dyear(year(Tpairs.ts_dstruct),      month(Tpairs.ts_dstruct),      day(Tpairs.ts_dstruct));

Tpairs = [Tpairs, table(tm_dyear,'VariableNames', {'tm_dyear'})];
Tpairs = [Tpairs, table(ts_dyear,'VariableNames', {'ts_dyear'})];
clear tm_dyear ts_dyear;



%% find unique
[Tpairs,ia,ic] = unique(Tpairs,'rows');

% find duplicates
id = setdiff(ia,ic);
if numel(id) > 0
    fprintf(1,'WARNING: found and deleted duplicate rows:\n');
    Tpairs(id,:);
end
Tpairs

% count
%npairs = numel(table2array(Tpairs(:,1)))
[nPairs, ndummy] = size(Tpairs);

%% select dates

% unpack Dates
for iPair=1:numel(datesFieldNames)
    datesFieldName1 = datesFieldNames{iPair}
    if isfield(DATES,datesFieldName1) == false
        switch iPair
            case 1  % origin epoch
                % round to nearest year
                yyyymmdd0 = nanmin([Tpairs.mast, Tpairs.slav]);
                yyyy=floor(yyyymmdd0/10000)*1000.
                yyyymmdd0 = yyyy + 0101;
                DATES.(datesFieldName1) = yyyymmdd0;
            case 2  % first acquisition date for stack
                yyyymmdd1 = nanmin([Tpairs.mast, Tpairs.slav]);
                DATES.(datesFieldName1) = yyyymmdd1;
            case 3  % last acquisition date for stack
                yyyymmdd2 = nanmax([Tpairs.mast, Tpairs.slav]);
                DATES.(datesFieldName1) = yyyymmdd2;
            otherwise
                error(sprintf('Unknown case i %d for date\n',iPair));
        end
    end
end
DATES


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


% %% work with logical variables
% isok1 = true(nPairs,1);
% isok1 = and(isok1,isfinite(Tpairs.mast));
% isok1 = and(isok1,isfinite(Tpairs.slav));
% isok1 = and(isok1,(Tpairs.mast >  DATES.yyyymmdd1));   % master date
% isok1 = and(isok1,(Tpairs.slav >  DATES.yyyymmdd1));   % slave  date
% isok1 = and(isok1,(Tpairs.mast <= DATES.yyyymmdd2));   % master date
% isok1 = and(isok1,(Tpairs.slav <= DATES.yyyymmdd2));   % slave  date
% % convert to indices
% iok1 = find(isok1 == true);

%% work with indices
iok = [1:nPairs];
iok = intersect(iok,find(isfinite(Tpairs.mast)));
iok = intersect(iok,find(isfinite(Tpairs.slav)));
iok = intersect(iok,find(Tpairs.mast >  DATES.yyyymmdd1));   % master date
iok = intersect(iok,find(Tpairs.slav >  DATES.yyyymmdd1));   % slave  date
iok = intersect(iok,find(Tpairs.mast <= DATES.yyyymmdd2));   % master date
iok = intersect(iok,find(Tpairs.slav <= DATES.yyyymmdd2));   % slave  date
Tpairs = Tpairs(iok,:);
[nPairs, ndummy] = size(Tpairs);


%% check that file exists and prune
kount = 0;
iok1 = zeros(nPairs,1);
Tpairs.filename
for iPair = 1:nPairs
    fname1 = char(Tpairs.filename{iPair});
    fprintf(1,'Grd file named fname1 = %s',fname1);
    if exist(fname1,'file') == 2
        kount = kount+1;
        iok1(iPair) = iPair;
        fprintf(1,' Exists.\n');
    else
        iok1(iPair) = 0;
        fprintf(1,' WARNING: Not found. Skipping.\n');
    end
end
iok=find(iok1 > 0);
Tpairs = Tpairs(iok,:);

% check size
[nPairs, ndummy] = size(Tpairs);
if nPairs ~= kount
    nPairs
    kount
    error(sprintf('miscount'));
end

%% Analyze graphs for several subsets of the data
%dataSets = {'ALL','MSF'};
%dataSets = {'ALL'};
%for idataSet = 1:numel(dataSets)
%
%dataSet1 = dataSets{idataSet}

if (nPairs > 0) == false
    error(sprintf('Too few pairs npairs = %d\n',nPairs));
else
    if CRITERIA.MSF
        nPasses = 2;
    else
        nPasses = 1;
    end
    for iPass = 1:nPasses
        
        %% get unique list of dates, sorted
        tu_cal   = unique([Tpairs.mast;  Tpairs.slav]);  % unique epochs
        tu_dyear = unique([Tpairs.tm_dyear, Tpairs.ts_dyear]);
        tu_dstruct = unique([Tpairs.tm_dstruct, Tpairs.ts_dstruct]);
        
        ntu = numel(tu_cal);
        tu_str = cell(ntu,1);
        for iPair=1:ntu
            tu_str{iPair} = sprintf('%8d',tu_cal(iPair));
        end
        
        
        %% list of edges by their calendar date YYYYMMDD
        Edges = [Tpairs.mast,Tpairs.slav];
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
            nPairs = npairs1;
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
                Weights = abs(Tpairs.bperp);
                ylab = 'Bperp';
                yunits = 'm';
            case 'Dt'  % useless because graph is a straight line
                % Use time span in days
                Weights = abs(Tpairs.dt);
                ylab = '\Delta t';
                yunits = 'days';
            case 'BpDt'
                %Use product of time span and Bperp
                %Weights = abs(zscore(Tpairs.bperp) .* zscore(Tpairs.dt));
                Weights = zscore(Tpairs.bperp) .* zscore(Tpairs.dt);
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
            
            %% find mininum spanning forest
            %if contains(dataSet, 'MSF') == true
            if CRITERIA.MSF && iPass == 1
                rs = abs(Weights);
                scalefactor = 1;
                [i_msf, i_rep, nCycles, QievMSF, TreesMSF, WeightsMSF] = findMininumSpanningForest(Qiev, itrees, Weights, scalefactor);
                nCycles
                % prune to MSF
                Tpairs = Tpairs(i_msf,:);
            end
            
            
            %% make a graph using Matlab functions
            DGRAPH = digraph(iedges(:,1),iedges(:,2),Weights,tu_str);
            Nodes = table2array(DGRAPH.Nodes);
            
            %% calculate importance of nodes
            importance = centrality(DGRAPH,'authorities','Importance',abs(DGRAPH.Edges.Weight));
            fprintf(1,'i,Nodes{i},tu_cal(i),importance(i) \n');
            for iPair=1:nepochs
                fprintf(1,'%3d %8s %8d %10.4f\n',iPair,Nodes{iPair},tu_cal(iPair),importance(iPair));
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
        
        
        [nPairs, ncolsmax] = size(Tpairs);
        
        %% check that arrays are OK
        % get information about first file
        %INFO=grdinfo3(char(grd_list(1)));
        % fname1 = TpairsOut.filename{1}
        % INFO1=grdinfo3(fname1);
        % [xgrd,ygrd,drho] = grdread3((fname1));
        % nf=nf+1;map_grd(fname1);
        
        
    end
    
    % if doPlots == 2
    %     pairListFilenameOutput = NaN;
    % else
    %     pairListFilenameOutput = sprintf('%s_%s_%3s_ipairs.txt' ...
    %         ,strrep(pairListFileNameInput,'.csv','') ...
    %         ,dataSet ...
    %         ,'select');
    %     writetable(Tpairs,pairListFilenameOutput);
    %     fprintf(1,'%s: successfully wrote new table named: %s\n',mfilename,pairListFilenameOutput);
    % end
    %
    % end
    
    return
end






