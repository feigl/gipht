function TSTP = generate_partials(DST,PST,TST)
% generate partials

fitfun = char(rowvec(PST.fitfun))
mparam = PST.mparam
ndata  = numel(DST.phaobs)
pnames = PST.names;

%copy initial values to final values
PST.p1 = PST.p0;

% initialize structure for speed
TSTP=struct('partial_wrt_1param',zeros(ndata,mparam));
fprintf(1,'ID, name, min, max, mean of partial derivative, delta (radians)\n');
for j=1:mparam
    %fprintf(1,'Calculating partial derivatives for %d %s\n',j,char(pnames{j}));
    % half a step down (left) in parameter
    PSTP1 = PST;
    PSTP1.p1(j) = PSTP1.p1(j) - 1.0d0 * PST.scale(j)/2.0;
    %PSTP1.p1(j) = -1.0d0 * PST0.scale(j);
    %ierr=check_struct(PSTP1);
    rng1 = feval(fitfun,DST,PSTP1,TST);
    %disp rng1; size(rng1)
    %fprintf(1,'Extrema for rng1 %g %g\n',min(rng1),max(rng1));
    % half a step up (right) in parameter
    PSTP2 = PST;
    PSTP2.p1(j) = PSTP2.p1(j) + 1.0d0 * PST.scale(j)/2.0;
    %PSTP2.p1(j) = +1.0d0 * PST0.scale(j);
    %ierr=check_struct(PSTP1);
    rng2 = feval(fitfun,DST,PSTP2,TST);
    %disp rng2; size(rng2)
    %fprintf(1,'%3d %s %g %g %g %g %g %g\n',j,char(pnames{j}),nanmin(rng1),nanmax(rng1),nanmin(rng2),nanmax(rng2),nanmin(rng2-rng1),nanmax(rng2-rng1));
    % partial derivative is difference (right minus left)
    der1 = colvec((rng2 - rng1) / PST.scale(j));
    %disp der1; size(der1)
    
    scl1 = PST.scale(j)*ones(ndata,1);
    der2 = zeros(ndata,1);
    deltap = ones(ndata,1);
    nbad = numel(find(isfinite(rng1)==0)) + numel(find(isfinite(rng2)==0));
    if nbad == 0
        for i=1:ndata
            derr=NaN;
            if abs(der1(i)) > 0.0
                % estimate partial derivate differently
                pp1 = PST.p1(j);
                fh = @(ppp)parameterfun18(ppp,i,j,fitfun,DST,PST,TST);
                [der2(i),derr, scl2(i)] = derivest(fh,pp1);
            else
                der2(i) = 0.0;
                scl2(i) = NaN;
            end
            
            if derr < 1.e-3
                ptlcol(i) = der2(i);
                deltap(i) = scl2(i);
            else
                ptlcol(i) = der1(i);
                deltap(i) = scl1(i);
            end
        end
    else
        ptlcol = zeros(ndata,1);
    end
    
    %     nbad = numel(find(isfinite(rng1)==0)) + numel(find(isfinite(rng2)==0));
    %     if nbad == 0
    %         ptlcol = der1;
    %     else
    %         ptlcol = zeros(ndata,1);
    %     end
    
    
    fprintf(1,'%3d %s %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n'...
        ,j,char(pnames{j})...
        ,nanmin(ptlcol),nanmax(ptlcol),nanmean(ptlcol)...
        ,nanmin(deltap),nanmax(deltap),nanmean(deltap));
    TSTP.partial_wrt_1param(1:ndata,j) = ptlcol;
end

return
end


