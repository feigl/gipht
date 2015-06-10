function status_message=write_table(filename,values,names,headers,ncols,iname)
%function status_message=write_table(filename,values,names,headers,ncols,iname)
% WRITE_TABLE: write an ascii table of data
% input parameters as column vectors
% return status
%
% ncols   = number of columns in table, including name
% header  = header lines to write
% iname   = column index for name (character string)
%           0 for no name
%           1 for first column
%           ncos for last column (no other columns allowed)
%
% Example:
%
% filename='test.dat'
% values=[1.0 2.0;1.1 2.1;0.9,1.9]
% names={'name1','name2','name3'}
% headers={'header 1','header 2'}
% ncols=3
% iname=3
% write_table =write_table(filename,values,names,headers,ncols,iname)
%
% write_table('test.dat',[1 2;3 4],{'name1','name2','name3'},{'header 1','header 2'},3,3)
%
% Will write a file named 'test.dat' with the following contents:
% * Written by feigl on 28-Sep-2011 at 09:22:17
% * header 1
% * header 2
%   +1.00000 +2.00000 name1
%   +1.10000 +2.10000 name2
%   +0.90000 +1.90000 name3
%
%
% Kurt Feigl, CNRS France April 96
% modified May 2002
% modified Feb 2003
% 2011-SEP-28 Kurt Feigl University of Wisconsin-Madison
% modified to use cell variables for names and headers
% 2012-MAR-01 fix bug with miscounted rows

if nargin ~= 6
%    disp 'Function: write an ASCII table with the right number of columns'
%    disp 'Usage: write_table(filename,values,names,headers,ncols,iname) '
   help(mfilename)
   error(sprintf('Wrong number of arguments. Need 6. Got %d',nargin));
end

status_message = sprintf('ERROR in write_table');

[fid,status_message] = fopen(filename,'w');
if fid == -1
  fprintf (1,'Error opening file %s\n',filename);
  disp (status_message)
  error
end

fmtstr = setstr(' ');

% how many rows?
%nrows = length(values(:,1));

% how many columns?
%ndim = length(values(1,:));

[nrows,ndim] = size(values);

% how many numerical values?
if iname == 0
   nvals = ncols;      % no names
else 
   nvals = ncols-1;    % names in first or last column
end

%whos names
if iscell(names) ~= 1
    warning('Array of names is not a cell. Converting');
    names0 = names;
    clear names;
    [nr,nc] = size(names0);
    nnames = max([nr,nc]);
    names = {};
    for i=1:nnames
        names{i} = names0(i,:);
    end    
    whos names
end

if (ndim ~= nvals) 
  warning(sprintf (1,'Mismatch in number of columns %d %d\n',ndim,nvals));
  nvals = ndim;
  %error (status_message)
  %return
end
if ismember(iname,[1,ncols]) == 1
    %if (nrows ~= numel(names))
    [ndummy,nnames] = size(names);
    nnames = max([ndummy,nnames]);
    %nnames = numel(names)
    %nnames = max(size(names{:}))
    if nrows ~= nnames
        warning(sprintf('Mismatch in number of rows %d %d\n',nrows,nnames));
        if nnames < nrows
            fprintf(1,'Padding.\n');
            for i=nnames+1:nrows
                names{i} = sprintf('row%05d',i);
            end
        end
%          whos names
%          %names{:}
%          names{1}
%          names{2}
    end
end


% make a format string
% names in first column
if iname == 1
   fmtstr=sprintf('%s',' %s');
end   
for i = 1:nvals
   % take standard deviation
   %sig = std(values(:,i));
   % mean
   %avg = mean(values(:,i));
   % extrema
   minv = min(abs(values(:,i)));
   maxv = max(abs(values(:,i)));
   sv = sort(values(:,i));
   dsv = sort(diff(sv));
   j = 1;
%   while dsv(j) <= 1.0d-99 & j < nrows-1
%     j = j+1;
%   end
   if j == nrows - 1 || maxv < 1.0d-15
      nfigs = 2;
      nleft = 1;
   else
%      mind = dsv(j);   
      if (abs(minv) > 0) 
         mind = minv;
      else
         mind = 1;
      end
      nfigs = ceil(log10(maxv/mind))+4;
      nleft = ceil(log10(maxv))+1;
   end


   if nleft >= 7
      ntot = nfigs + 5;
      ndec = nfigs;
      lett = 'E';
   elseif nleft >= 3
      ntot = nfigs + 4;
      ndec = nfigs - nleft + 3;
      lett = 'f';
   elseif nleft <= -2
      ntot = nfigs + 5;
      ndec = nfigs;
      lett = 'E';
   else
      ntot = nfigs + 3;
      nleft = 0;
      ndec = nfigs;
      lett = 'f';
   end

   if ~(ntot > 0 & ntot < 20) | ~(ndec >= 0 & ndec < 20)
      ntot = 12;
      ndec =  4;
      lett = 'e';
   end
   
   fmtstr=sprintf('%s%s%d.%d%s',fmtstr,' %+0#',ntot,ndec,lett);
end
% names in last column
if iname == ncols
   fmtstr=sprintf('%s%s',fmtstr,' %s');
end   

% how many headers?
% nheaders = length(headers);
% [nheaders,n] = size(headers);
nheaders = numel(headers);

% write header lines
if nheaders > 0
    s = date;
    c = clock;
    fprintf(fid,'* Written by %s on %s at %2.2d:%2.2d:%02.0f\n'...
        ,getenv('USER'),s,c(4),c(5),c(6));   
    for i = 1:nheaders
        %fprintf(fid,'* %s\n',headers(i,:));
        fprintf(fid,'* %s\n',char(headers{i}));
        ferror(fid);
    end
end

for i = 1:nrows
    switch iname
        % no names
        case 0
            [s,status_message] = sprintf(fmtstr,values(i,1:nvals));
            
            % names in first column
        case 1
            %     [s,status_message] = sprintf(fmtstr,names(i,:),values(i,1:nvals));
            [s,status_message] = sprintf(fmtstr,char(names{i}),values(i,1:nvals));
            
            % names in last column
        case ncols
            %[s,status_message] = sprintf(fmtstr,values(i,1:nvals),names(i,:));
            %fprintf(1,'names{%d} is %s\n',i,char(names{i}));
            %size(char(names{i}))
%            [s,status_message] = sprintf(fmtstr,values(i,1:nvals),char(names{i}));
            [s,status_message] = sprintf(fmtstr,values(i,1:nvals),names{i});
        otherwise
            disp(status_message);
            error(status_message);
    end
    if length(status_message)
        disp(status_message);
        error(status_message);
    end
    fprintf (fid,'%s\n',s);
    ferror(fid);
end

status_message = sprintf ('Wrote %d lines to %s\n',i+nheaders+1,filename);
fclose(fid);


return;
