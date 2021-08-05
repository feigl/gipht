function AR = build_AR(A,d,a,c,X,Y,rowstep)
%% build A'*R where R is data spatial covariance matrix (too large to store)
% inputs 
% A == design matrix
% d == data vector
% X == X coordinates of pixels in data vector [meters]
% Y == Y coordinates of pixels in data vector [meters]
% true representation
% references:
% From: Kurt Feigl <feigl@wisc.edu>
% Date: Tuesday, December 15, 2020 at 2:44 PM
% To: "Reinisch, Elena Cristina" <ecreinisch@lanl.gov>
% Subject: Re: [EXTERNAL] Matlab code for Geostatistical inversion
% Would you please confirm that I am looking at the correct piece of code (attached)?
% https://raw.githubusercontent.com/feigl/gipht/Brady/HTC_mesh_build/mesh_range2d_HTC.m

% From: Reinisch, Elena Cristina <ecreinisch@lanl.gov>
% Sent: Tuesday, December 15, 2020 18:03
% To: Kurt Feigl <feigl@wisc.edu>
% Subject: Re: [EXTERNAL] Matlab code for Geostatistical inversion 
% That file is correct!


%c = 0.01*(0.005/dt)^2; %3e-6; %5e-8;
%a = 230; %m
% step_size = 3000
[arrows, arcols1] = size(A');
arcols = numel(d);

if arcols1 ~= arcols
    arcols1
    arcols
    error(sprintf('Number of columns is incorrect.\n'));
else
    AR = zeros(size(A'));
    fprintf(1,'Starting %s...\n',mfilename);
    start_time = tic;
    
    % row index
    rind = [(0:rowstep:arcols), arcols]
    
    %% original code is a bit cryptic
    % rind = [(0:3000:arcols), arcols];
    % for i = 1:numel(rind)-1
    %     AR(:, rind(i)+1:rind(i+1)) = A'*[c*exp(-3.*(abs(sqrt((repmat(datx_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(datx_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2+(repmat(daty_vec(rind(i)+1:rind(i+1))', [arcols, 1]) - repmat(daty_vec(:), [1, numel(rind(i)+1:rind(i+1))])).^2))./a))];
    % end
    %         ARfac = 1;
    %         AR = AR*ARfac;
    
    % this should be equivalent
    nind = numel(rind)-1;
    for i = 1:nind
        X1 = repmat( X(rind(i)+1:rind(i+1))', [arcols, 1]);
        X2 = repmat( X(:)                   , [1, numel(rind(i)+1:rind(i+1))]);
        Y1 = repmat( Y(rind(i)+1:rind(i+1))', [arcols, 1]);
        Y2 = repmat( Y(:)                   , [1, numel(rind(i)+1:rind(i+1))]);
        T1 = exp( -3.*(abs(sqrt( (X1 - X2).^2 + (Y1 - Y2).^2)/a)));
        P1 = c*(A'*T1);
        %     whos
        %     disp('c');size(c)
        %     disp('At');size(A')
        %     disp('T1');size(T1)
        %     disp('P1');size(P1)
        %     disp('size(AR(:, rind(i)+1:rind(i+1)))');size(AR(:, rind(i)+1:rind(i+1)))
        AR(:, rind(i)+1:rind(i+1)) = P1;
        if mod(i,10) == 1
            fprintf(1,'i = %6d of %6d Elapsed time in seconds: %.1f\n',i,nind,toc(start_time));
        end
    end
end


return
end

