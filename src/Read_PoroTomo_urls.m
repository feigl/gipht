

clear all; close all; clc

% list of files available on GDR
%File_url{1} = 'https://gdr.openei.org/files/828/Brady_Obs_Metadata.csv';Short_name{1}='VIBE';
% File_url{2} = 'https://gdr.openei.org/files/827/Reftek_metadata.csv';
% File_url{3} = 'https://gdr.openei.org/files/826/Nodal_continuous_metadata%20(1).csv';

%File_names{1} = {'Brady_Obs_Metadata.csv'};
%File_names{2} = {'Reftek_metadata.csv'};
%File_names{3} = {'Nodal_continuous_metadata%20(1).csv'};

% for i = 1:numel(File_url)
%     urlwrite(File_url{i},char(File_names{i}));
% end

% these files are local versions after removing ^M characters by hand AND
% removing commas by Hand.
File_names{1} = '/data/PoroTomo/Metadata/Brady_Obs_MetadataEOL.csv';        Short_name{1}='VIBE';
File_names{2} = '/data/PoroTomo/Metadata/Reftek_metadataEOL.csv';           Short_name{2}='REFT';
File_names{3} = '/data/PoroTomo/Metadata/Nodal_continuous_metadataEOL.csv'; Short_name{3}='NODE';

% build a structure
META = struct

for i = 1:numel(File_names)
    META.(Short_name{i}) = csv2struct(char(File_names{i}))
    %META.(Short_name{i}) = csv2struct2(char(File_names{i}))

end


    






