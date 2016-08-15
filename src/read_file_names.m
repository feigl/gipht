function [pfnames, mdate, imast, sdate, islav, hamb, ddays, t1, t2, idatatypes, mpercys] = read_file_names(file_name)
% read the phase file names and decimal years for master and slave
% FILE_NAMES.DAT format
%  decimal.year_master  decimal.year_slave  phase_file_name  % comments
% 1993.5589  1995.4438  psp_10575_20438_ort.pha  % 1993 JUL 24 1995 JUN 12 
% 1993.5589  1995.7315  psp_10575_21941_ort.pha  % 1993 JUL 24 1995 SEP 25
% [pfnames, mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_file_names(file_name)
% [1xNcel, 1xNcel, Nx1dbl, 1xNcel, Nx1dbl, Nx1dbl, Nx1dbl, Nx1dbl, Nx1dbl] where N=#pairs
% for use with DIAPASON DTOOLS
% 20160524 handle file names from GMT5SAR
fprintf(1,'%s begins ...\n',mfilename);

%% check to see that file exists
if fexist(file_name) == 0 
    error(sprintf('Cannot open file named %s containing names of phase files.\n',file_name));
    return
end

%% decide to read GMT grid files or not
if numel(strfind(file_name,'grd')) > 0
    file_type = 2;
else
    file_type = 1;
end

%% open file
fid = fopen(file_name,'r');
if fid <= 0 
    error(sprintf('Cannot open file named %s containing names of phase files.\n',file_name));
    return
end


switch file_type
    case 1
        %% old format with decimal years
        % 1992.59836065574 1993.46027397260 ../intf/1992220_1993169/phasefilt_ll.grd     % wrapped phase in radians
        % 1992.59836065574 1993.46027397260 ../intf/1992220_1993169/unwrap_mask_ll.grd   % unwrapped phase in radians
        
        % Read the arguments into a Cell Array
        CDAT=textscan(fid,'%f%f%s','CommentStyle','%');
        fclose(fid);
        [nrows,ncols] = size(CDAT);
        if ncols == 3 && nrows > 0
            t1    = CDAT{1};  % First column of dat file is the Master Date
            t2    = CDAT{2};  % Second column of dat file is the Slave Date
            names = CDAT{3};  % Third column of dat file is the phase file names
            np = nrows;
            idatatypes = zeros(np,1);
        else
            error(sprintf('miscount with file_type = %d',file_type))
        end
    case 2
        %% new format with GMT grid files
        %  ../intf/1992220_1993169/phasefilt_ll.grd     % wrapped phase in radians
        %  ../intf/1992220_1993169/unwrap_mask_ll.grd   % unwrapped phase in radians
        
        % Read the arguments into a Cell Array
        CDAT=textscan(fid,'%s %s %d %f %s','CommentStyle','%');
        fclose(fid);
        [nrows,ncols] = size(CDAT);
        if ncols == 5 && nrows > 0
            yyyymmdd1    = CDAT{1};  % Master Date
            yyyymmdd2    = CDAT{2};  % Slave Date
            idatatypes   = CDAT{3};  % data type
            mpercys      = CDAT{4};  % fringe spacing meters per cycle
            names        = CDAT{5};  % third column of dat file is the grid file name
            
            np = nrows;
            for i=1:np
                str1 = char(yyyymmdd1(i));
                str2 = char(yyyymmdd2(i));
                y1 = str2num(str1(1:4));
                y2 = str2num(str2(1:4));
                m1 = str2num(str1(5:6));
                m2 = str2num(str2(5:6));
                d1 = str2num(str1(7:8));
                d2 = str2num(str2(7:8));
                
                
                %[y1,m1,d1] = yeardoy2yyyyymmdd(yyyy1,doy1);
                %[y2,m2,d2] = yeardoy2yyyyymmdd(yyyy2,doy2);
                % t1 = datetime(y1,m1,d1,'TimeZone','UTC');t1.Format = 'yyyy-MM-dd_HH:mm:SSSSSS [ZZZZ]'  ; % 24 hour clock
                % t2 = datetime(y2,m2,d2,'TimeZone','UTC');t2.Format = 'yyyy-MM-dd_HH:mm:SSSSSS [ZZZZ]'  ; % 24 hour clock
                t1(i) = datetime(y1,m1,d1,'TimeZone','UTC');t1.Format = 'yyyy-MM-dd'  ; % 24 hour clock
                t2(i) = datetime(y2,m2,d2,'TimeZone','UTC');t2.Format = 'yyyy-MM-dd'  ; % 24 hour clock
                                
%                 % try to guess data type from file name
%                 if numel(strfind(str1,'phasefilt')) > 0
%                     idatatypes(i) = 0;
%                 elseif numel(strfind(str1,'qgradx')) > 0
%                     idatatypes(i) = -1;
%                 elseif numel(strfind(str1,'unwrap_mask')) > 0
%                     idatatypes(i) = 2;
%                 else
%                     errror(sprintf('cannot parse data type from file name %s\n',str1));
%                 end
            end
        else
            error(sprintf('miscount with file_type = %d',file_type))
        end
    otherwise
        error(sprintf('Unknown case: file_type is %d\n',file_type));
end

%% make column vectors
idatatypes = colvec(idatatypes)
mpercys = colvec(mpercys)

% % Check if any pairs were read
% np = numel(n0);
% if np < 1
%     fprintf(1,'ERROR No lines in %s \n',file_name);
%     error(sprintf('Cannot read phase file names dat file named %s',file_name)); 
%     return
% end
% 
% % Check that three args were read for each pair
% if np ~= numel(t0) || np ~= numel(t1)
%     fprintf(1,'ERROR formatting problem in %s \n',file_name);
%     error 'Check phase file names file format' 
%     return
% end
fprintf(1,'Read %4d pairs of phase files \n',np);

% Now building dummy arguments to match interferogram.lst code
uhamb = zeros(1,np)';
uimast = zeros(1,np)';
uislav = zeros(1,np)';
uddays = zeros(1,np)';
for N = 1:np
	ipair = N;
	if N== 1
		umdate = cellstr(sprintf('PAIR%03d-MAST',N));
		usdate = cellstr(sprintf('PAIR%03d-SLAV',N));
	else
		umdate = {umdate{:} sprintf('PAIR%03d-MAST',N)};
		usdate = {usdate{:} sprintf('PAIR%03d-SLAV',N)};
	end
	uhamb(N) = ipair*10;
% 	uimast(N) = -1000 * ipair;
% 	uislav(N) = +1000 * ipair;
	uimast(N) = +1000 * ipair + 1;
	uislav(N) = +1000 * ipair + 2;
	uddays(N) = 100 * ipair;
end
% pfnames = n0';
pfnames = names;
mdate = umdate;
sdate = usdate;
hamb = uhamb;
imast = uimast;
islav = uislav;
ddays = uddays;

return
