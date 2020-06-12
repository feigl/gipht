function nreset = set_comsol_parameters2(model,PST)
% transfer PST.p1 of parameters from PST structure to MPH file
% 20200611 Kurt Feigl

narginchk(2,2);

mparams = numel(PST.p1);
verbose = 1;
nreset = 0;

fprintf(1,'Entering %s\n',mfilename);
PST.p1

model = mphload(PST.fileNameMPH);

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
            if verbose == 1
                fprintf(1,'    Setting %d %-8s to %12.4e %s\n',ip,PST.name{ip},PST.p1(ip),PST.description{ip});
            end
            %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
            model.param.set(sprintf('%s',PST.name{ip}),sprintf('%e[%s]',PST.p1(ip),PST.units{ip}),PST.description{ip});
            nreset = nreset + 1;
        end
        %     end
        %
        %
        %         % % set Poisson's ratio
        %         % %PST.NUP=0.2;
        %         % if isfinite(PST.NUP)==1
        %         %     model.param.set('NUP', sprintf('%e',PST.NUP),'Poissons Ratio');
        %         % end
        %         %
        %         % % set Young's Modulus
        %         % if isfinite(PST.EYM)==1
        %         %     model.param.set('EYM', sprintf('%e[Pa]',PST.EYM), 'Youngs Modulus');
        %         % end
        %         %
        %         % model.param.set('LENG', sprintf('%e[m]',PST.LENG), 'ellipsoid a axis');
        %         % model.param.set('WIDTH', sprintf('%e[m]',PST.WIDTH), 'ellipsoid b axis');
        %         % %model.param.set('PRES', sprintf('%e[Pa]',PST.PRES), 'pressure');
        %
        %         %end
        %     end
        % end
        
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
        mphsave(model,PST.fileNameMPH,'optimize','speed');
    %
    %     if verbose == 1
    %         fprintf(1,'Finished task in %.1f seconds\n',toc(tstart));
    %     end
    %

return
end



