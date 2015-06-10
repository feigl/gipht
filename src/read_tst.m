function TST = read_tst(fnametst)
%function TST = read_tst(fnametst)
% read a TST structure from a flat ASCII file named fnametst
% partial derivative of range (in cycles) 
% with respect to parameter (in units specified by scale factor PST.scale)

% open file for data
fid = fopen(fnametst,'r');
if fid <= 0
    error(sprintf('Cannot open TST file called %s\n',fnametst));
end

% number of columns
%ncols = 19;
[nsize, nread] = fscanf(fid,'%d %d\n',2);
nrows = nsize(1)
ncols = nsize(2)

fprintf(1,'Number of columns (paramaters) in file is %d\n',ncols);

% count the number of lines
nlines = -1;
aline = ' ';
while aline ~= -1 
    %[A, nread] = fscanf(fid,fmt,ncols);
    aline = fgetl(fid);
    if aline ~= -1
       nlines = nlines+1;
    end
end
fprintf(1,'Number of lines (data points) in file is %d\n',nlines);
frewind(fid);


% write a format statement
fmt = '%g'; % 
% must match number of real variables in fprintf statement below
for i=1:ncols-1
    fmt = sprintf('%s %s',fmt,'%g');
end
fmt = sprintf('%s\\n',fmt);

% initialize structure for speed
TST.partial_wrt_1param        = zeros(nlines,1); % partial derivative

aline = fgetl(fid); % skip header line
k = 0;
nread = ncols;
aline = ' ';
while aline ~= -1
    %[A, nread] = fscanf(fid,fmt,ncols);
    aline = fgetl(fid);
    if aline ~= -1
        [A, nread] = sscanf(aline, fmt, ncols);
        if nread == ncols
            tstart = tic;
            k = k+1;
            TST.partial_wrt_1param(k,1:ncols)            = rowvec(A(1:ncols)); % partial
        end
         if mod(k,1000) == 0
            fprintf(1,'Read line %5d %12.4g\n',k,toc(tstart));
        end
    end
end
fclose(fid);
TST

return

