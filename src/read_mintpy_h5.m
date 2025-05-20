function [INFO, ATTR, DATA] = read_mintpy_h5(h5_file_name)
%% Read data from HDF5 file in Matlab

% Ref https://www.mathworks.com/help/matlab/high-level-functions.html

% 2021/06/11 Kurt Feigl
% 2021/10/18 Kurt Feigl
% 2024/09/09 Kurt Feigl
% 2025/04/14 Finally get this correct

narginchk(1,1);
fprintf(1,'reading solution from MintPy in HDF5 file named %s\n', h5_file_name);

%% Read attributes from HDF5 file
ATTR=struct;
DATA=struct;
INFO = h5info(h5_file_name)
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
    if numel(strfind(name1,'X_STEP'))>0 || numel(strfind(name1,'Y_STEP')) >0
        fprintf('%s is %f\n',char(name1),ATTR.(name1));
    end
end

ndatasets=numel(INFO.Datasets)

%% read data sets
for i=1:ndatasets
    name1 = INFO.Datasets(i).Name;    
    fprintf(1,'Looking for dataset named %s .... \n',name1)
    DATA.names{i}=name1;
end
DATA.names

for i=1:ndatasets
    name1 = INFO.Datasets(i).Name;    
    fprintf(1,'Starting to read dataset %s .... \n',name1)

    A = h5read(h5_file_name, sprintf('//%s',name1));
    %whos A
    ndim = numel(size(A));
    fprintf(1,'Array %s has dimension = %d numel = %d\n',name1,ndim,numel(A));
    if ndim == 1
            [nr_dat,nc_dat] = size(A);
            fprintf(1,'nrows = %d by ncols = %d\n',nr_dat,nc_dat);
            DATA.(name1) = colvec(A);

            if strcmp(name1,'date')
                % handle dates as Matlab datetime
                DATA.datetime=datetime(DATA.(name1),'Format','yyyyMMdd','TimeZone','UTC');
            end
    elseif ndim == 2
        %[nr_dat,nc_dat] = size(A);
        [nrA,ncA] = size(A);
        fprintf(1,'from file nrows = %d by ncols = %d\n',nrA,ncA);


        if nrA == 1 || ncA == 1
            % vector
            DATA.(name1) = colvec(A);
            if strcmp(name1,'date')
                % handle dates as Matlab datetime
                DATA.datetime=datetime(DATA.(name1),'Format','yyyyMMdd','TimeZone','UTC');
            end
        else 
            % handle 2-D array
            % 2025/04/14
            DATA.(name1) = transpose(A);
            [nr_dat,nc_dat] = size(DATA.(name1));
            fprintf(1,'after transposing nrows = %d by ncols = %d\n',nr_dat,nc_dat);

            % make arrays of coordinates
            if  isfield(ATTR,'Y_STEP') && isfield(ATTR,'Y_FIRST') ...
                    && isfield(ATTR,'X_STEP') && isfield(ATTR,'X_FIRST')
                DATA.vecX=[ATTR.X_FIRST : ATTR.X_STEP: ATTR.X_FIRST+(ATTR.WIDTH  -1)*ATTR.X_STEP]';
                fprintf(1,'number of X coordinates = %d\n',numel(DATA.vecX));
                DATA.vecY=[ATTR.Y_FIRST : ATTR.Y_STEP: ATTR.Y_FIRST+(ATTR.LENGTH -1)*ATTR.Y_STEP]';
                fprintf(1,'number of Y coordinates = %d\n',numel(DATA.vecY));
                % 2025/04/14
                [DATA.XGRD,DATA.YGRD]=meshgrid(DATA.vecX,DATA.vecY);
                [nr_grd,nc_grd] = size(DATA.XGRD);
                if nr_grd ~= nr_dat
                    nr_dat
                    nr_grd
                    error('number of rows do not agree');
                elseif  nc_grd ~= nc_dat
                    nc_dat
                    nc_grd
                    error('number of columns do not agree');
                end
            else
                warning('in data set named %s of cannot find coordinates of 2-D grid',name1, h5_file_name);
            end
        end
    elseif ndim == 3
            [nrA3, ncA3, n_layers] = size(A);
            fprintf(1,'nrows = %d by ncols = %d by nlayers = %d\n',nrA3,ncA3,n_layers);
            B = nan(ncA3,nrA3,n_layers);
            for j=1:n_layers
                A1 = squeeze(A(:,:,j));
                % 2025/04/14
                B(:,:,j)=transpose(A1);
            end
            DATA.(name1) = B;
    else
            error(sprintf('Cannot handle dimension %d\n',ndim));
    end   
end
%% update names of data sets
dataSetNames=fieldnames(DATA);
return
end



