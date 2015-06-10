function [pfnames, mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_file_names(file_name)
% read the phase file names and decimal years for master and slave
% FILE_NAMES.DAT format
%  decimal.year_master  decimal.year_slave  phase_file_name  % comments
% 1993.5589  1995.4438  psp_10575_20438_ort.pha  % 1993 JUL 24 1995 JUN 12 
% 1993.5589  1995.7315  psp_10575_21941_ort.pha  % 1993 JUL 24 1995 SEP 25
% [pfnames, mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_file_names(file_name)
% [1xNcel, 1xNcel, Nx1dbl, 1xNcel, Nx1dbl, Nx1dbl, Nx1dbl, Nx1dbl, Nx1dbl] where N=#pairs
% for use with DIAPASON DTOOLS
fprintf(1,'%s begins ...\n',mfilename);

fid = fopen(file_name,'r');
if fid <= 0 
    error 'Cannot open phase file_names descriptor'
    return
end

% Read the arguments into a Cell Array
	idat=textscan(fid,'%f%f%s','CommentStyle','%');
fclose(fid);

t0 = idat{1};     % First column of dat file is the Master Date
t1 = idat{2};     % Second column of dat file is the Slave Date
n0 = idat{3};     % Third column of dat file is the phase file names

% Check if any pairs were read
np = numel(n0);
if np < 1
    fprintf(1,'ERROR No lines in %s \n',file_name);
    error(sprintf('Cannot read phase file names dat file named %s',file_name)); 
    return
end

% Check that three args were read for each pair
if np ~= numel(t0) || np ~= numel(t1)
    fprintf(1,'ERROR formatting problem in %s \n',file_name);
    error 'Check phase file names file format' 
    return
end
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
pfnames = n0;
mdate = umdate;
sdate = usdate;
hamb = uhamb;
imast = uimast;
islav = uislav;
ddays = uddays;

return
