function ierr = write_gin(PST,fnameout)
%function ierr = write_gin(PST,fnameout)
% Write a PST structure containing parameters into .gin file with name fname

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
    % 1st header line 
    fprintf(fid,'Parameter_Name Initial_Value PlusMinusBound\n');

    % write one line per parameter
    fmt = '%32s %20.10E %20.10E\n';     
    kount = 0;
    for i=1:PST.mparam
        % Do not write NaNs or derived parameters
        %&& (strcmp(char(PST.flag(i),'E#') > 0 || strcmp(char(PST.flag(i)),'F#') > 0)) ...
        %if numel(PST.names(i)) > 2 ...           
        if      isfinite(PST.p1(i)) == 1 ...
                && isfinite(PST.sigma(i)) == 1
            fprintf(fid,fmt,char(PST.names{i}),PST.p1(i),PST.sigma(i));
            kount = kount+1;
        end
    end
    % check count
    %if kount == PST.mparam
    if kount > 0 && kount <= PST.mparam
        fprintf(1,'Wrote %d fields for %d parameters to file %s\n',nfn,kount,fnameout);
        fclose(fid);
        ierr = 0;
    else
        ierr = 1;
         warning(sprintf('mismatch between kount (%d) and PST.mparam (%d)\n',kount,PST.mparam));
    end
end
return
end
