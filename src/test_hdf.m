%% set up HDF5 - not needed with Matlab release 2024B
% download pre-made binaries from https://www.hdfgroup.org/download-hdf5/

% 2026/03/27 -- new macbook pro with Apple Silicon chip
% '/Applications/HDF_Group/HDF5/2.1.0/lib/plugin');

% sudo mkdir /usr/local/hdf5
% cd /usr/local/hdf5/
% ls -l /Applications/HDF_Group/HDF5/2.1.0
%sudo ln -s /Applications/HDF_Group/HDF5/2.1.0/* .

setenv('HDF5_PLUGIN_PATH','/Applications/HDF_Group/HDF5/2.1.0/lib/plugin');
%setenv('HDF5_PLUGIN_PATH','/usr/local/hdf5/lib/plugin')

% give up on above
% https://formulae.brew.sh/formula/hdf5
% brew hdf5
hdf5_plugin_path = getenv('HDF5_PLUGIN_PATH')
dir(hdf5_plugin_path)
%cd(hdf5_plugin_path)


%% get information
for i=1:2
    switch i
        case 1 % no compression  -- works
            fname='/Volumes/feigl/insar/BRADY/SDK/DESCENDING/mintpy63/avgPhaseVelocity.h5'
            data_set = '/velocity'
            % mamba activate mintpy
            % info.py /Volumes/feigl/insar/BRADY/SDK/DESCENDING/mintpy63/avgPhaseVelocity.h5
            % HDF5 dataset "/velocity            ": shape=(114, 90)           , dtype=float32   , compression=None
        case 2 % lzf compression -- fails
            fname='/Volumes/feigl/insar/BRADY/SDK/DESCENDING/mintpy63/inputs/geometryGeo.h5'
            data_set = '/height'
            % info.py /Volumes/feigl/insar/BRADY/SDK/DESCENDING/mintpy63/inputs/geometryGeo.h5
            % HDF5 dataset "/height              ": shape=(114, 90)           , dtype=float32   , compression=lzf
            % HDF5 dataset "/incidenceAngle      ": shape=(114, 90)           , dtype=float32   , compression=lzf
            % HDF5 dataset "/slantRangeDistance  ": shape=(114, 90)           , dtype=float32   , compression=lzf
            % HDF5 dataset "/waterMask           ": shape=(114, 90)           , dtype=bool      , compression=lzf
        otherwise
            error('unknown case %d',i)
    end
   
    I=h5info(fname)
    I.Datasets
    I.Datasets.Name
    I.Attributes
    nDataSets=numel(I.Datasets)
    for j=1:nDataSets
        fprintf(1,'%d %s\n',j,I.Datasets(j).Name);
    end
    %% read data
    D = h5read(fname,data_set);
    whos D
end

%
% Error using h5readc
% The HDF5 library encountered an error and produced the following stack trace information:
%
%
%     H5PL_load    can't find plugin. Check either HDF5_VOL_CONNECTOR, HDF5_PLUGIN_PATH, default location, or path set
%     by H5PLxxx functions
% Error in h5read (line 95)
%     [data,var_class] = h5readc(Filename,Dataset,start,count,stride);
%                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Error in test_hdf (line 38)
% D = h5read(fname,'/height')


return

