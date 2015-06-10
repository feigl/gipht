function ierr = print_parameters_nicely(PST,fname)
%function ierr = print_parameters_nicely(PST)
% print a table of parameters with non-zero adjustments 

ierr = 0;

if exist('fname','var') == 1
    fd = fopen(fname,'a');
else
    fd = 0;
end


% unload parameters
mparam = PST.mparam; % number of parameters
pnames = PST.names;  % parameter names
p0 = PST.p0;         % initial estimate
p1 = PST.p1;         % final estimate
%pflags = PST.flag;   % flag
psig   = PST.sigma;  % 1-sigma uncertainty
upb    = PST.ub;     % upper bound
lpb    = PST.lb;     % lower bound

imode = 1;           % short format string

% print a nice header
nwide = 58;
aline='_';
for i=1:nwide-1
    aline=strcat(aline,'_');
end
fprintf(1,'%s\n',aline);
hedfmt = '%-32s  %-10s   %-10s\n';
fprintf(1,hedfmt,'Parameter','Estimate','Uncertainty');
fprintf(1,'%s\n',aline);
if fd > 0
    fprintf(fd,'%s\n',aline);
    fprintf(fd,hedfmt,'Parameter','Estimate','Uncertainty');
    fprintf(fd,'%s\n',aline);
end

for j=1:mparam
    adj = p1(j)-p0(j);
    if abs(adj) > 0
        outfmt = getfmt(p1(j),pnames{j},imode);
        
        fprintf(1      ,outfmt,pnames{j} ,p1(j),psig(j));
        if fd > 0
            fprintf(fd ,outfmt,pnames{j} ,p1(j),psig(j));
        end
    end
end
fprintf(1,'%s\n',aline);
if fd > 0
    fprintf(fd,'%s\n',aline);
    fclose(fd);
end

return


