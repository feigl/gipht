function phamod = parameterfun18(ppp,irow,jcol,fitfun,DST,PST,TST)
% anonymous function call to go with DERIVEST

    %fprintf(1,'In %s\n',mfilename);
    
    
    % extract one line (a record for one observation point) from DST
    DST = extract_dst(DST,irow);
    
    % check that DST structure is OK
    ierr = check_struct(DST);
    
    % call fitting function first time to initialize
    % temporary storage structure TST    
    [rng,TST] = feval(fitfun,DST,PST);
    
%     disp ppp; whos dsp
%     disp ppp; size(ppp)
%     phamod = zeros(size(ppp));
    for i=1:numel(ppp)
        % reset the variable parameter to be differentiated
        PST.p1(jcol) = ppp(i);
        
        % call the function
        phamod(i) = feval(fitfun,DST,PST,TST);
    end
    return
end

