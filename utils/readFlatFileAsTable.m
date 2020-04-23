function Tflat = readFlatFileAsTable(FileName,nHeaderLines)
%function Tflat = readFlatFileAsTable(FileName,nHeaderLines)
% read flat text file as a Matlab table
% inputs:
%    FileName    == name of input file
%    HeaderLines == number of header lines to skip [default = 1]
% outputs:
%    Tflat       == Table with variable names
% 20200402 Kurt Feigl

if nargin < 2
    nHeaderLines = 1;
end

% get first header line
fid = fopen(FileName,'rt');
headerline = fgetl(fid);
fclose(fid);
% change commas to single spaces
headerline = strrep(headerline,',',' ');
varnames = matlab.lang.makeValidName(split(headerline));
varnames = strrep(varnames,'___','_');
varnames = strrep(varnames,'__','_');
varnames = strrep(varnames,'x_','');
nvars = numel(varnames);


% read table
Tflat = readtable(FileName, 'FileType','text','ReadVariableNames', 0 ...
    ,'HeaderLines',nHeaderLines ...
    ,'MultipleDelimsAsOne',1);

% if delimiter is not specified, then guess
% , 'Delimiter', {'space','tab'},
%     ,'Delimiter', 'space' ...

% count rows and columns
[nrows, ncols] = size(Tflat);

% truncate if necessary
if nvars > ncols
    fprintf(1,'Warning: file contains more variables (nvars = %d) than columns (ncols = %d). Ignoring last %d columns.\n'...
        ,nvars,ncols,nvars-ncols);
    varnames = varnames(1:ncols);
end
% pad if necessary
if nvars < ncols
    fprintf(1,'Warning: file contains fewer variables (nvars = %d) than columns (ncols = %d). Padding last %d columns.\n'...
        ,nvars,ncols,ncols-nvars);
    for i=nvars+1:ncols
        varnames{i} = sprintf('Variable%03d',i);
    end
end

Tflat.Properties.VariableNames = varnames;
fprintf(1,'Read Table %s with %d rows and %d columns and following variable names:\n',FileName,nrows,ncols);
for i=1:numel(varnames)
    fprintf(1,'%s ',char(varnames{i}));
end
fprintf(1,'\n');

return
end




