function ierr = write_pst(PST,fnameout)
%function ierr = write_pst(PST,fnameout)
% Write a PST structure containing parameters with name fname

if nargin ~= 2
    error('Missing datafilename');
end

nfn = numel(fieldnames(PST));
if nfn == 10
    fprintf(1,'Found %d fields in PST. Setting scale factor to 1.\n',nfn);
    for i=1:PST.mparam
        PST.scale(i) = 1.0;
    end
%elseif nfn == 11
elseif nfn == 12
    %fprintf(1,'Found %d fields in PST structure.\n',nfn);
else
    error(sprintf('Wrong number of fields in PST. Found %d . Expected 12\n',nfn));
end

fid = fopen(fnameout,'w');
if fid <= 0
    error(sprintf('Cannot open data file called %s\n',fnameout));
    ierr = 1;
else
    % 1st header line contains number of lines and the name of the file
    fprintf(fid,'%d %s %s\n',PST.mparam,char(PST.fitfun),char(PST.datafilename));
    
    
    % 2nd header line contains the names of the columns
%     varnames = fieldnames(PST);
%     for i=1:numel(varnames);
%         fprintf(fid,' %s',varnames{i});
%     end
%     fprintf(fid,'\n');
    fprintf(fid,'PST.i(i),char(PST.names{i}),PST.p0(i),PST.p1(i),PST.lb(i),PST.ub(i),PST.sigma(i),PST.scale(i)\n');

    
    % write a format statement
    ncols = 8; % must match number of variables in fprintf statement below
    fmt = '%4d %s';
    for i=1:ncols-2
        fmt = sprintf('%s %s',fmt,'%20.10e');
    end
    fmt = sprintf('%s\\n',fmt);
    
    kount = 0;
    for i=1:PST.mparam
%        fprintf(fid,'%f',PST.i,char(PST.names{i}),PST.p0(i),PST.p1(i),PST.lb(i),PST.ub(i),PST.sigma(i));
%        fprintf(fid,fmt,PST.i(i),char(PST.names{i}),PST.p0(i),PST.p1(i),PST.lb(i),PST.ub(i),PST.sigma(i));
        fprintf(fid,fmt,PST.i(i),char(PST.names{i}),PST.p0(i),PST.p1(i),PST.lb(i),PST.ub(i),PST.sigma(i),PST.scale(i));
        kount = kount+1;
    end
    % check count
    if kount == PST.mparam
        fprintf(1,'Wrote %d fields for %d parameters to file %s\n',nfn,kount,fnameout);
        fclose(fid);
        ierr = 0;
    else
        ierr = 1;
        error(sprintf('kount (%d) not equal to PST.mparam (%d)\n',kount,PST.mparam));
    end
end
return
end
