function [rng0,TSTP] = generate_partials(DST,PST,TST)
% generate partials

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

fitfun = char(rowvec(PST.fitfun));
mparam = PST.mparam;
ndata  = numel(DST.phaobs);
pnames = PST.names;

%copy initial values to final values
PST.p1 = PST.p0;

% evaluate exact fitting function at initial estimate of parameters
rng0 = feval(fitfun,DST,PST,TST);
DST.phamod = rng0;


% Jacobian matrix is partial derivative of observable with respect to parameter
dDdP = zeros(ndata,mparam);
%fprintf(1,'ID, name, min, max, mean of partial derivative, delta (radians)\n');

%ptlcol = zeros(ndata,1);

% count free parameters by considering difference between bounds
db = PST.ub-PST.lb;
ifree = find(abs(db) > 0.);
mfree = numel(ifree);

% % decide if we need to run comsol
% do_comsol = 0;
% if numel(strfind(lower(char(PST.datafilename)),'.mph')) > 0
%     % data file exists
%     if fexist(char(PST.datafilename)) == 1
%         do_comsol = 1;
%     end
% end



%parfor j=1:mparam % does not work if calling comsol
for j=1:mparam
    % set all elements to zero
    ptlcol = zeros(ndata,1);

    if db(j) > 0.
    fprintf(1,'--- Calculating partial derivatives for parameter %5d (of %5d) %s %12.4e\n',j,mparam,char(pnames{j}),PST.scale(j));
    
%     if do_comsol == 1
%         model(j) = mphload(PST.datafilename);
%         mphnew = sprintf(strrep(char(PST.datafilename),'.mph',sprintf('_%03d.mph',j)));
%         mphsave(model(j),mphnew);
%         addAttachedFiles(gcp,model(j));
%     end
    
    %     Error using generate_partials>(parfor supply) (line 41)
    %     An UndefinedFunction error was thrown on the workers for 'model'.  This might be because the file containing 'model' is not accessible on the
    %     workers.  Use addAttachedFiles(pool, files) to specify the required files to be attached.  See the documentation for
    %     'parallel.Pool/addAttachedFiles' for more details.
    %         Error in generate_partials (line 41)
    %         parfor j=1:mparam % does not work if calling comsol
    %             Error in test_partials (line 13)
    %             [rngdum,TSTP] = generate_partials(DST,PST,TST);
    %             Caused by:
    %             Undefined function or variable "model".
    
    % half a step down (left) in parameter
    PSTP1 = PST;
    PSTP1.p1(j) = PSTP1.p1(j) - 1.0d0 * PST.scale(j)/2.0;
    PSTP1.flag{j} = 'F#';
%     if do_comsol == 1
%         PSTP1.datafilename = mphnew;
%     end
    rng1 = feval(fitfun,DST,PSTP1,TST);
    
    % half a step up (right) in parameter
    PSTP2 = PST;
    PSTP2.p1(j) = PSTP2.p1(j) + 1.0d0 * PST.scale(j)/2.0;
    PSTP2.flag{j} = 'F#';
%     if do_comsol == 1
%         PSTP2.datafilename = mphnew;
%     end
    rng2 = feval(fitfun,DST,PSTP2,TST);
    
    % partial derivative is difference (right minus left)
    der1 = colvec((rng2 - rng1) / PST.scale(j));
    %der1 = colvec((rng2 - rng1) / db / 2.0);
    
    % find valid elements
%     iok1 = find(isfinite(der1)==1);  % finite value
%     iok2 = find(abs(der1)>0.0);      % non-zero
%     iok  = intersect(iok1,iok2);     % both of above
    
    iok = 1:ndata;
    
    
    % overwrite with valid elements
    ptlcol(iok) = der1(iok);
    
%     figure;hold on;
%     set(gcf,'DefaultTextInterpreter','none');
%     plot(rng1,'b');
%     plot(rng2,'r');
%     plot(der1,'k-');
%     xlabel('index');
%     ylabel('partial derivative');
%     legend('rng1','rng2','der1');
%     title(sprintf('parameter %d %s\n',j,char(pnames{j})));
    
    % count bad elements
    nbad = ndata - numel(iok);
    if nbad ~= 0
        warning(sprintf('Replaced %d partial derivatives with zero.\n',nbad));
    end
    
    
    %     fprintf(1,' %03d %s %12.5E %12.5E %12.5E %12.5E\n'...
    %         ,j,char(pnames{j})...
    %         ,nanmin(ptlcol),nanmax(ptlcol),nanmean(ptlcol));
    %
    
    end
    % return partial derivative wrt 1 parameter as 1 column in TSTP structure
    dDdP(:,j) = ptlcol;
end

fprintf(1,'\nFinished generating partial derivatives for %d free parameters\n',mfree);

% % display non-zero elements of the Hessian matrix of partial derivatives
% figure;
% spy(TSTP.partial_wrt_1param);
% xlabel(sprintf('mparam = %d columns',mparam));
% ylabel(sprintf('ndata = %d rows',ndata));

TSTP.partial_wrt_1param = sparse(dDdP);


fprintf(1,'\n\n----------------   %s ends   at %s ----------\n',upper(mfilename),datestr(now,31));

return
end


