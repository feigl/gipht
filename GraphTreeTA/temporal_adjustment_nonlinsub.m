function [mse] = temporal_adjustment_nonlinsub(data,tm,ts,tbreaks,tfunc,metaparams, x, V)
% sub function of temporal adjustment for quick inversion and estimation of
% mse.  Used in nonlinear optimization (fmincon)
% 
% Elena C. Baluyut; UW-Madison
% 2015-11-09

%% Initialize 
% Find number of data and define variables
ndat = numel(data);
tu = colvec(sort(unique([tm; ts]))); % unique set of epochs
me = numel(tu);  % number of epochs
metaparams(find(metaparams == inf)) = x;

%% Build design matrix G
        mparams = numel(time_function(tfunc, min([tm; tu]), tbreaks, metaparams));
        
        % Set up design matrix G
        G = zeros(ndat,mparams);
        % Build design matrix G from temporal function
        for i=1:ndat
            tfm = time_function(tfunc, tm(i), tbreaks, metaparams); % time function for master
            tfs = time_function(tfunc, ts(i), tbreaks, metaparams); % time function for slave
            for jj=1:mparams
                G(i,jj) = tfs(jj)-tfm(jj); % time function for pair
            end
        end % loop over data


    

 % All other parameterizations     
      [pest, psig, mse, Vx,] = ls_with_cov(G, data, V);
%   end
return



