%% read an MPH file output by COMSOL
% 20200420 Kurt Feigl

%% initialize
close all
clear variables
verbose = 0
if verbose == 1
    echo on;
    diary(sprintf('%s_%s.log',mfilename,datestr(now,30)));
else
    echo off;
end
    

%% Comsol server must be started
try
    mphstart
catch ME
    ME.identifier
    warning('COMSOL mph server is not started.\nTrying to restart');
    switch computer
        case 'MACI64'
           [status, output] = system('/Applications/COMSOL54/Multiphysics/bin/comsol mphserver &')
        otherwise
            error('Unknown computer');
    end
end

%% import some utilities
import com.comsol.model.util.*

%% check on status
status = which('mphload')    
if contains(status,'mphload.p') == 1
    fprintf(1,'Matlab Livelink server is available.\n');

    fprintf(1,'If this is new to you, consider the following commands:\n');
    fprintf(1,'  help mphload\n');
    fprintf(1,'  mphnavigator\n');
    fprintf(1,'  mphgetproperties\n');
    fprintf(1,'  mphdoc(\''mphload\'')\n');
else
    error('Matlab Livelink server is not available.\n');
end

% get a list of the appropriate MPH files
% 01/LdM_3DFSI_spheroid_Tvisco_OPT2019DEC02WORKS_01.mpho'
 
DIRstruct=dir('*/*.mpho');

% loop over MPH files 
nFiles = numel(DIRstruct)
for i=1:nFiles        
    fileNameMPH = strcat(DIRstruct(i).folder,filesep,DIRstruct(i).name)
    
    
    % load the MPH file
    model=mphload(fileNameMPH);
    
    % get some information
    info0 = mphsolutioninfo(model);
    
    % make the available plots
    MSTRUCT = mphmodel(model,'-STRUCT');
    for j=1:numel(MSTRUCT.result)
        plotGroup1 = MSTRUCT.result{j}
        if numel(plotGroup1) > 0
            figure;
            try
                mphplot(model,plotGroup1);
                fileNamePDF = strrep(fileNameMPH,'.mpho',sprintf('_%s.pdf',plotGroup1));
                title(plotGroup1);
                print(gcf,'-dpdf',fileNamePDF,'-r600','-fillpage','-painters');
            catch ME
                ME
            end
        end
    end
    
    % get a table
    % we have to know its tag!
    Tstruct = mphtable(model,'tbl4');
    % get the data
    ArrayData = Tstruct.data;
    [nrows,ncols] = size(ArrayData);
    % get the headers and make them into legitimate variable names
    varnames = Tstruct.headers;
    varnames = matlab.lang.makeValidName(varnames);
    % make a table
    Ttable = array2table(ArrayData);   
    % assign names to variables
    Ttable.Properties.VariableNames = varnames;
    
    % Add variables to table
    Ttable = [Ttable, table(i*ones(nrows,1),'VariableNames', {'ISolution'})];

    % summarize table
    summary(Ttable);

    % write the table for this run
    fileNameCSV = strrep(fileNameMPH,'.mpho','.csv')
    writetable(Ttable,fileNameCSV); 
    
    % build big table
    if i==1
        TBIG = Ttable;
        ncols1 = ncols;
        nrows1 = nrows;
    else
        if ncols == ncols1
            TBIG(nrows1+1:nrows1+nrows,:) = Ttable;
            nrows1 = nrows1+nrows;
        else
            nrows
            nrows1
            ncols
            ncols1
            error('miscount');
        end
    end
end
fileNameBIG = sprintf('%s_BIG.csv',mfilename)
writetable(TBIG,fileNameBIG);


