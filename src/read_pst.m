function PST = read_pst(fnamein,datafilename)
%function PST = read_pst(fnamein)
% read structure PST for parameters from file named fnamein

fid = fopen(fnamein,'r');
if fid <= 0
    error(sprintf('Cannot open data file called %s\n',fnamein));
    ierr = 1;
else
    % read the first line
    C = textscan(fid,'%s%s%s\n',1);
    mparam = str2double(C{1});
    fitfun = C{2};
    datafilename = C{3};
    
    % enter these into the structure
    PST.mparam       = mparam;
    PST.datafilename = datafilename;
    PST.fitfun       = fitfun;
   
    % read the second line, containing the names of the columns
    % ignore them for now
    C = textscan(fid,'%s\n',1);


    % write a format statement
    ncols = 9; % must match number of variables in fprintf statement below
    fmt = '%4d%s';
    for i=1:ncols-2
        fmt = sprintf('%s%s',fmt,'%f');
    end
    fmt = sprintf('%s\\n',fmt);
    for i = 1:mparam
        C = textscan(fid,fmt,1);
        if numel(C) ~= ncols
            error(sprintf('Number of parameters listed in file (%d) does not equal number specified (%d)\n',numel(C),ncols));
        end
        j = 1;
        PST.i(i)      = C{j};j=j+1; % index of parameter
        PST.names(i)  = C{j};j=j+1; % name of parameter
        PST.p0(i)     = C{j};j=j+1; % initial value of parameter
        PST.p1(i)     = C{j};j=j+1; % current best, or final value of parameter
        PST.lb(i)     = C{j};j=j+1; % lower bound
        PST.ub(i)     = C{j};j=j+1; % upper bound
        PST.sigma(i)  = C{j};j=j+1; % uncertainty
        PST.scale(i)  = C{j};j=j+1; % scale factor
        PST.flag(i)   = C{j};j=j+1; % flag
    end
    
    fn=fieldnames(PST);
    % make
    for i=1:numel(fn)
        F1=getfield(PST,fn{i});
        [nr,nc] = size(F1);
        % make row vectors into column vectors
        if nr == 1
            F1=reshape(F1,nc,1);
            PST=setfield(PST,fn{i},F1);
        end
        [nr,nc] = size(F1);
        %fprintf(1,'DST %s is %d rows by %d columns.\n',fn{i},nr,nc);
    end

    %     kount = 0;
    %     for i=1:mparams
    %
    %         fprintf(fid,fmt,i,pnames{i},p(i),bounds(i,1),bounds(i,2));
    %         kount = kount+1;
    %     end
    return
    fclose(fid);
end