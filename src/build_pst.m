function PST = build_pst(fitfun,mparam,p0,p1,psig,pnames,bounds,datafilename,pscl,pflag)
%function ierr = write_fitfunin(fnamein,fitfun,mparams,p,pnames,bounds,pscl)
% Build a PST structure containing parameters with name fname
% Kurt 20101111 add scale factor 


if nargin ~= 10
    error('Wrong number of arguments');
end

% enter these into the structure
PST.mparam       = mparam;
PST.datafilename = datafilename;
PST.fitfun       = fitfun;

for i = 1:mparam
    PST.i(i)      = i;           % index of parameter
    PST.names{i}  = pnames{i};   % name of parameter
    PST.p0(i)     = p0(i);       % initial value of parameter
    PST.p1(i)     = p1(i);       % current best, or final value of parameter
    PST.lb(i)     = bounds(i,1); % lower bound
    PST.ub(i)     = bounds(i,2); % upper bound
    PST.sigma(i)  = psig(i);     % uncertainty   
    PST.scale(i)  = pscl(i);     % scale factor in same units as above
%    PST.flag{i}   = 'N#';        % flag
    PST.flag{i}   = pflag{i};     % flag
end

fn=fieldnames(PST);

% make column vectors, except for names
%for i=1:numel(fn)
for i=[1 3:numel(fn)]
    F1=getfield(PST,fn{i});
    [nr,nc] = size(F1);
    % make row vectors into column vectors
    if nr == 1
        F1=reshape(F1,nc,1);
        PST=setfield(PST,fn{i},F1);
    end
    [nr,nc] = size(F1);
    %fprintf(1,'PST %s is %d rows by %d columns.\n',fn{i},nr,nc);
end
nfn = numel(fieldnames(PST));
%fprintf(1,'Wrote %d fields in PST\n',nfn);

return
end