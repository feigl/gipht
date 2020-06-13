function nreset = set_comsol_parameters2(PST)
% transfer PST.p1 of parameters from PST structure to MPH file
% 20200611 Kurt Feigl

narginchk(1,1);

%fprintf(1,'Entering %s\n',mfilename);

mparams = numel(PST.p1);
verbose = 1;
nreset = 0;

% extract model from structure
model = PST.model;

% for i=1:mparams
%     fprintf(1,'%10.4g ',PST.p1(i));
% end
% fprintf(1,'\n');

%model = mphload(PST.fileNameMPH);

for ip=1:mparams  % loop over indices to comsol parameters
    if isfinite(PST.p1(ip)) == 1
        %         if verbose == 1
        %             fprintf(1,'%d %10s %12.4e %s\n',ic,pnamesCS{ic},PST.p1(ic),descrs{ic});
        %         end
        % get the index of the old value
        %it = find(strcmp(Tparams.name,PST.name(ip)) == true);
        %pold = Tparams.value(it);
        pold = nan;
        
        %if strfind(pnamesCS{ic},'PRES') > 0
        % get index to gipht parameter
        %ip = get_parameter_index('CS_PRES_in_Pa_pressure__________',PST.names);
        %ip = get_parameter_index(sprintf('CS_%s*',pnamesCS{ic}),PST.names);
        %p1 = PST.p1(ip); % initial estimate of gipht parameter
        db = PST.ub(ip)-PST.lb(ip); % range of bounds in gipht parameter
        
        
        if  db > eps
            %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
            pname1 = sprintf('%s',PST.name{ip});
            unit1 = char(model.param.evaluateUnit(pname1));
            desc1 = char(model.param.descr(pname1));
            if verbose == 1
                %fprintf(1,'    Setting %d %-8s to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),PST.description{ip});
                fprintf(1,'    Setting %d %-8s to %12.4e [%s] %s\n',ip,pname1,PST.p1(ip),unit1,desc1);
            end
            model.param.set(pname1,sprintf('%e[%s]',PST.p1(ip),unit1),desc1);
            nreset = nreset + 1;
        end        
    end
end   
    %     if verbose == 1
    %         nreset
    %         fprintf(1,'Saving %s\n',mphname);
    %     end
    %     tstart = tic;
    %
    %     %% Save the parameters to the MPH file
    %     %mphsave(model,sprintf('%s.mph',mphname));
    %     %mphsave(model,mphname);
    %
    %     % 20140829 need to instance it
    %     %model.save(mphname);
    %
        %20150901 call function
    %   mphsave(model,char(mphname));
        % 20200612 call function, optimizing for speed
%        mphsave(model,PST.fileNameMPH,'optimize','speed');
    %
    %     if verbose == 1
    %         fprintf(1,'Finished task in %.1f seconds\n',toc(tstart));
    %     end
    %

return
end



