function nreset = set_comsol_parameters(PST, mphname, pnamesCS, values, dims, descrs)
% transfer values of parameters from PST to MPH file

verbose = 0;

if verbose == 1;
    fprintf(1,'Loading %s\n',mphname);
end
tstart=tic;
model = mphload(mphname);
if verbose == 1
    fprintf(1,'Finished task in %.1f seconds\n',toc(tstart));
end

nreset = 0;
for ic=1:numel(values)  % loop over indices to comsol parameters
    if isfinite(values(ic)) == 1
%         if verbose == 1
%             fprintf(1,'%d %10s %12.4e %s\n',ic,pnamesCS{ic},values(ic),descrs{ic});
%         end
        
        
        %if strfind(pnamesCS{ic},'PRES') > 0
        % get index to gipht parameter
        %ip = get_parameter_index('CS_PRES_in_Pa_pressure__________',PST.names);
        ip = get_parameter_index(sprintf('CS_%s*',pnamesCS{ic}),PST.names);
        p1 = PST.p1(ip); % initial estimate of gipht parameter
        db = PST.ub(ip)-PST.lb(ip); % range of bounds in gipht parameter
        
        
        if (abs(p1) > 0)
            if abs(PST.p1(ip)-values(ic)) > eps*abs(values(ic))
                %if verbose == 1
                    fprintf(1,'    Re-setting %d %-8s %-32s from %12.4e to %12.4e %s\n',ip,pnamesCS{ic},PST.names{ip},values(ic),PST.p1(ip),descrs{ic});
                %end
                %model.param.set('PRES', sprintf('%e[Pa]',PST.p1(ip)),descrs{ic});
                model.param.set(sprintf('%s',pnamesCS{ic}),sprintf('%e%s',PST.p1(ip),dims{ic}),descrs{ic});
                nreset = nreset + 1;
            end
        end
        
        
        % % set Poisson's ratio
        % %PST.NUP=0.2;
        % if isfinite(PST.NUP)==1
        %     model.param.set('NUP', sprintf('%e',PST.NUP),'Poissons Ratio');
        % end
        %
        % % set Young's Modulus
        % if isfinite(PST.EYM)==1
        %     model.param.set('EYM', sprintf('%e[Pa]',PST.EYM), 'Youngs Modulus');
        % end
        %
        % model.param.set('LENG', sprintf('%e[m]',PST.LENG), 'ellipsoid a axis');
        % model.param.set('WIDTH', sprintf('%e[m]',PST.WIDTH), 'ellipsoid b axis');
        % %model.param.set('PRES', sprintf('%e[Pa]',PST.PRES), 'pressure');
        
        %end
    end
end



if verbose == 1;
    nreset
    fprintf(1,'Saving %s\n',mphname);
end
tstart = tic;
%mphsave(model,sprintf('%s.mph',mphname));
%mphsave(model,mphname);

% 20140829 need to instance it 
model.save(mphname);

if verbose == 1
    fprintf(1,'Finished task in %.1f seconds\n',toc(tstart));
end

return




