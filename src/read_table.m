function [values, name] = read_table(filename,nheaders,ncols,iname)
%function [values,name]=read_table(filename,nheaders,ncols,iname)
% read an ascii table of data
% return parameters as column vectors
%
% ncols   = number of columns in table, including name
% nheader = number of header lines to skip
% iname   = column index for name (character string)
%           0 for no name
%           1 for first column
%           ncols for last column (no other columns allowed)
%
% Kurt Feigl, CNRS France, June     2004
%
% WARNING! This version works for MATLAB version 6.1
% earlier versions of MATLAB may require earlier versions of this routine!
% Updated for Matlab R2102b Kurt Feigl 2013/10/01

values = NaN;
name = sprintf('___');

if nargin < 1
    disp 'Function: read an ASCII table with Matlab v. 6.1'
    disp ' '
    disp 'Usage: read_table(filename,nheader,ncols,iname) '
    disp ' '
    disp 'nheader = number of header lines to skip'
    disp 'ncols   = number of columns in table, including name'
    disp 'iname   = column index for name (character string)'
    disp '          0 for no name'
    disp '          1 for first column'
    disp '          ncols for last column (no other columns allowed)'
    disp ' '
    disp 'Example: [xy,names]=read_table(''fig5.xyr'',2,6,6)'
    disp ' '
    error('Wrong number of arguments.');
end

if nargin < 2
    nheaders = 0;
    warning('Assuming no headers!');
end

if nargin < 3
    ncols = 2;
    warning('Assuming 2 columns!');
end

if nargin < 4
    iname = 0;
    warning('Assuming no names!');
end



fid = fopen(filename,'r');

fmtstr = sprintf(' ');

% how many numerical values?
if iname == 0
    nvals = ncols;      % no names
else
    nvals = ncols-1;    % names in first or last column
end
values = zeros(1,nvals);

% make a format string
for i = 1:nvals
    fmtstr=sprintf('%s %s','%f',fmtstr);
end

% skip header lines
for i = 1:nheaders
    %   [s,count] = fscanf(fid, '%[^\n]s',1);
    s = fgetl(fid);
    fprintf(1,'%s\n',s);
end

in = 0;
count = 1;
while count > 0
    % get one line
    s = fgetl(fid);
    
    % check that line is a string
    if isstr(s)
        in = in+1;
        
        switch iname
            case 0
                name1 = sprintf('Line %10d',in);
            case 1  % skip over name in first column
                s=strjust(s,'left');
                iblanks=strfind(s,' ');
                iblank1=min(iblanks);
                j = max([2,iblank1]);
                name1 = s(1:j-1);
                %s=s(j:length(j));
                % 20131001
                s = s(j:end);
            case ncols  % get ascii names in last column
                if length(s) > 1
                    name1 = lststr(s);
                else
                    name1 = '_';
                end
            otherwise
                error(sprintf('Unknown iname = %d\n',iname));
        end
        
        % load character array
        name(in,1:length(name1)) = name1;
        
        % read numerical items
        [dumarr,count,errmsg] = sscanf(s,fmtstr,nvals);
        if errmsg
            disp(errmsg)
        end;
        if (count == nvals)
            values(in,1:nvals) = dumarr';
        end;
    else
        count = 0;
    end;
end;

% 
fprintf(1,'\nNumber lines read: %d\n',in);
fclose(fid);

return;
end




function s=lststr(t)

% find the last string s starting with a blank in t


m = length(t);

% last item is a white space
if isspace(t(m))
    while isspace(t(m)),
        m = m -1;
    end
end;

while ~ isspace(t(m)),
    m = m -1;
end

s = t(m:length(t));

return;
end
