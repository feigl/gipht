function kselect = select_pixels(pselect,npix,grdd,mskd)
%function kselect = select_pixels(pselect,npix,grdd,mskd)
% choose pixels according to pselect, return 1-D vector in list

switch pselect
    case 0
        fprintf(1,'Selecting all finite pixels.\n');
        kkeep = find(isfinite(grdd)==1);
    case {1,3} % Choose pixels randomly
        fprintf(1,'Selecting %d pixels randomly\n',npix);
        %20130624 - fprintf(1,'Initializing random number generator.\n');
        %        Replace the default stream with a stream whose seed is based on CLOCK, so
        %        RAND will return different values in different MATLAB sessions.  NOTE: It
        %        is usually not desirable to do this more than once per MATLAB session.
        % 20130624 - Do something about the error message below.
        %         RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
        %         Warning: The RandStream.setDefaultStream static method will be removed in a future
        %         release.  Use RandStream.setGlobalStream instead.
        %         > In RandStream.RandStream>RandStream.setDefaultStream at 456
        %         In gipht_step1 at 691
        %         In gipht at 119
        %          randi Pseudorandom integers from a uniform discrete distribution.
        %     R = randi(IMAX,N) returns an N-by-N matrix containing pseudorandom
        %     integer values drawn from the discrete uniform distribution on 1:IMAX.
        %     randi(IMAX,M,N) or randi(IMAX,[M,N]) returns an M-by-N matrix.

        % 20170830 WRONG: kkeep =  randi(npix,[numel(grdd),1]); 
        % 20170830 RIGHT: draw npix values from interval [1,nx*ny]
        kkeep =  randi([1,numel(grdd)],npix,1);
        fprintf(1,'Saving %d pixel indices to kkeep.mat \n',npix);
        save kkeep kkeep;
    case 2      %20130624 - My understanding is that if we do NOT initialize, then the RNG
        %will generate a different number in each session of GIPHT.
        % See if file exists
        fdk = fopen ('kkeep.mat','r'); 
        if fdk == -1 
            error(sprintf('Could not open kkeep.mat .\nPlease set select = 1.\n'));
        else
            fclose(fdk);
            fprintf(1,'Reading previous pixel indices from existing file kkeep.mat.\n');
            load kkeep; 
            if numel(kkeep) ~= npix 
                error(sprintf(...
                    'Number of pixels (%d) does not match number in kkeep.mat.\nPlease set select = 1.\n',npix));
            end
        end
    case {5,7}
        % mask 
        kkeep=find(isfinite(grdd .* mskd)==1);      
    case {4,6,9}
        warning(sprintf('pselect == %d not implemented\n',pselect));
    otherwise,
        error(sprintf('Unkown value of pselect (%d)\n',pselect));
end

whos

%% select finite data 
kfinite = find(isfinite(grdd)==1);

%% compute intersection
kselect = intersect(kkeep,kfinite);

return;

