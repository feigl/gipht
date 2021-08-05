function [INFO, ATTR, DATA] = read_mintpy_h5(h5_file_name)
%% Read data from HDF5 file in Matlab

% Ref https://www.mathworks.com/help/matlab/high-level-functions.html

% 2021/06/11 Kurt Feigl

%% Read attributes from HDF5 file
ATTR=struct;
INFO = h5info(h5_file_name);
[num_atr, ~] = size(INFO.Attributes);
for i = 1:num_atr
    fprintf(1,'%s %s\n',char(INFO.Attributes(i).Name), char(INFO.Attributes(i).Value));
    name1 = INFO.Attributes(i).Name;
    name1=  matlab.lang.makeValidName(name1);
    valu1 = INFO.Attributes(i).Value;
    % try decoding string
    valu2 = str2double(valu1);   
    if isfinite(valu2)
        ATTR.(name1) = valu2;  % value is numeric
    else
        ATTR.(name1) = valu1;  % value is string
    end   
end

% az_lks      = h5readatt(vel_file,'/','ALOOKS');             az_lks      = str2double(az_lks{1});
% rg_lks      = h5readatt(vel_file,'/','RLOOKS');             rg_lks      = str2double(rg_lks{1});
% head_angle  = h5readatt(vel_file,'/','HEADING');            head_angle  = str2double(head_angle{1});
% asc_desc    = h5readatt(vel_file,'/','ORBIT_DIRECTION');
% prf         = h5readatt(vel_file,'/','PRF');                prf         = str2double(prf{1});
% az_pix_size = h5readatt(vel_file,'/','AZIMUTH_PIXEL_SIZE'); az_pix_size = str2double(az_pix_size{1});
% rg_pix_size = h5readatt(vel_file,'/','RANGE_PIXEL_SIZE');   rg_pix_size = str2double(rg_pix_size{1});
% ref_lat     = h5readatt(vel_file,'/','REF_LAT');            ref_lat     = str2double(ref_lat{1});
% ref_lon     = h5readatt(vel_file,'/','REF_LON');            ref_lon     = str2double(ref_lon{1});
% wvl         = h5readatt(vel_file,'/','WAVELENGTH');         wvl         = str2double(wvl{1});
% lon1        = h5readatt(vel_file,'/','X_FIRST');            x0          = str2double(x0{1});
% lat1        = h5readatt(vel_file,'/','Y_FIRST');            y0          = str2double(y0{1});
% lon_step    = h5readatt(vel_file,'/','X_STEP');             x_step      = str2double(x_step{1});
% lat_step    = h5readatt(vel_file,'/','Y_STEP');             y_step      = str2double(y_step{1});

lon_vec=[ATTR.X_FIRST : ATTR.X_STEP : ATTR.X_FIRST+ATTR.X_STEP*ATTR.WIDTH];
lat_vec=[ATTR.Y_FIRST : ATTR.Y_STEP : ATTR.Y_FIRST+ATTR.Y_STEP*ATTR.LENGTH];



ndatasets=numel(INFO.Datasets);

%% read data sets
for i=1:ndatasets
    name1 = INFO.Datasets(i).Name
    dataSetNames{i} = name1;
    A = h5read(h5_file_name, sprintf('/%s',name1));
    ndim = numel(size(A));
    fprintf(1,'Array %s has dimension = %d numel = %d ',name1,ndim,numel(A)); 
    switch ndim
        case 1
            [nrows,ncols] = size(A);
             fprintf(1,'nrows = %d by ncols = %d\n',nrows,ncols);
        case 2
            [nrows,ncols] = size(A);
             fprintf(1,'nrows = %d by ncols = %d\n',nrows,ncols);
        case 3
            [nrows, ncols, nlayers] = size(A);
            fprintf(1,'nrows = %d by ncols = %d by nlayers = %d\n',nrows,ncols,nlayers);
        otherwise
            error(sprintf('Cannot handle dimension %d\n',ndim));
    end
    
    if ATTR.Y_STEP < 0
        % transpose and flip up-to-down
        switch ndim
            case 1
                DATA.(name1) = colvec(A);
            case 2
                DATA.(name1) = flipud(transpose(A));
            case 3
                B = nan(ncols,nrows,nlayers);
                for j=1:nlayers
                    A1 = squeeze(A(:,:,j));
                    B(:,:,j)=flipud(transpose(A1));
                end
                DATA.(name1) = B;
                otherwise
            error(sprintf('Cannot handle dimension %d\n',ndim));
        end
    else
        % transpose only
         switch ndim
            case 1
                DATA.(name1) = colvec(A);
            case 2
                DATA.(name1) = transpose(A);
            case 3
                B = nan(ncols,nrows,nlayers);
                for j=1:nlayers
                    A1 = squeeze(A(:,:,j));
                    B(:,:,j)=transpose(A1);
                end
                DATA.(name1) = B;
            otherwise
                error(sprintf('Cannot handle dimension %d\n',ndim));
        end
        
    end
end
if ATTR.Y_STEP < 0
    lat_vec = lat_vec(end:-1:1);
end

% make grid of coordinates

jcols = repmat([1:ATTR.WIDTH]  ,  ATTR.LENGTH,          1);
irows = repmat([1:ATTR.LENGTH]',            1, ATTR.WIDTH);

fieldNames=fieldnames(DATA);
if sum(contains(fieldNames,'longitude'))
   DATA.longitudes  = lon_vec(jcols);  % do not overwrite existing field
else
   DATA.longitude = lon_vec(jcols);
end
if sum(contains(fieldNames,'latitude'))
   DATA.latitudes  = lat_vec(irows);  % do not overwrite existing field
else
   DATA.latitude = lat_vec(irows);
end

DATA.lon_vec = lon_vec; % row vector
DATA.lat_vec = transpose(lat_vec); % column vector

% update names of data sets
dataSetNames=fieldnames(DATA);
return
end



