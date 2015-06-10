function csv_filename = convert_abaqus_dat_to_csv(dat_filename)
%function csv_filename = convert_abaqus_dat_to_csv(dat_filename)
% calls stepstripper routine to parse output from ABAQUS .DAT file 
% writes a Comma Separated Values (CSV) file

% [pathstr,fname,ext] = fileparts(dat_filename)
% csv_filename = fullfile(pathstr,filesep,fname,'csv')
% tbl_filename = fullfile(pathstr,filesep,fname,'tbl')

tbl_filename = strrep(dat_filename,'.dat','.tbl');
csv_filename = strrep(dat_filename,'.dat','.csv');


commandline = sprintf('%s %s %s %s\n'...
    ,get_executable_name('stepstripper.cpp')...
    ,dat_filename ...
    ,csv_filename ...
    ,tbl_filename ...
    );
[status, result] = system(commandline);

if status ~= 0
    status
    result
    error('call to system failed');
end


if fexist(csv_filename) ~= 1
    error(sprintf('Could not find CSV file named %s\n',csv_filename));
end

return;


