function var = read_insar_data(isat,itrack,varname)
%function var = read_insar_data(satname,itrack,varname)
% given satellite track id, return east component of unit vector pointing from target on ground to sensor in orbit along line of sight
% suitable for use with COMSOL
% example:
%      satname = 'ALOS'
%      itrack = 112
%      varname = 'RangeChangeInMeters'
%      var = read_insar_data(satname,itrack,varname)
% 20191108 Kurt Feigl
switch satname
    case 1 % 'ALOS'
        switch itrack
            case 112
                %filename = 'MAU2_TS_AT112.csv';
                filename = 'MAU2_TS_AT112.xlsx';
            otherwise
                error(sprintf('Unknown itrack %d\n',itrack));
        end
    otherwise
        error(sprintf('Unknown satname %d\n',satname));
end                

T=readtable(filename);
S=table2struct(T,'ToScalar',true)
if isfield(S,varname)
    var = S.(varname);
else
    var = NaN;
    error(sprintf('Requested variable %s is not a column in data file named %s\n',varname,filename));
end
return
end