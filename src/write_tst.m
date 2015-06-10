function ierr = write_tst(TST,fitfun,fnametst)
%function ierr = write_tst(TST,fitfun,fnametst)
% write the partial derivative matrix to a file

% open file for TST
fid = fopen(fnametst,'w');
if fid <= 0
    error(sprintf('Cannot open TST file called %s\n',fnametst));
    return
end

% number of columns
%ncols = 20; 
varnames = fieldnames(TST);
%ncols = numel(varnames);

% number of rows
%ndata = numel(TST.partial_wrt_1param)
[ndata,ncols] = size(TST.partial_wrt_1param)
TST

fprintf(fid,'%d %d %s\n',ndata,ncols,fitfun);


% write a format statement
fmt = '%20.10E'; % 
% must match number of real variables in fprintf statement below
for i=1:ncols-1
    fmt = sprintf('%s %s',fmt,'%20.10E');
end
fmt = sprintf('%s\\n',fmt)
numel(fmt)
size(fmt)
size(fmt)/8

% reserve space for the model
phamod = 0.;
% loop over rows
kount = 0;
for i=1:ndata
    kount = kount+1;
    fprintf(fid,fmt,TST.partial_wrt_1param(i,1:ncols));
end

% check count
if kount == ndata
    fprintf(1,'Wrote meta data for %d pixels to file %s\n',kount,fnametst);
    fclose(fid);
    ierr = 0;
else
    error(sprintf('kount (%d) not equal to ndata (%d)\n',kount,ndata));
    ierr = 1;
end

return

