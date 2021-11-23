function LIMITS = get_site_dims(sitecode5)
%function LIMITS = get_site_dims(sitecode5)
% return bounding box 
% 2021/10/18 Kurt Feigl
home=getenv('HOME');
fname=strcat(home,filesep,'site_dims.txt');
fid=fopen(fname,'rt');
if fid > 0
    while 1
        tline = fgetl(fid);
        if ischar(tline)
            % look for 6 characters with a colon
            if numel(tline) == 6 && numel(strfind(tline,':')) > 0
                % compare first 5 characters
                if strncmpi(tline,sitecode5,5)
                    for i=1:3
                        tline=fgetl(fid);
                        tline=strrep(tline,'-R','');
                        tline=strrep(tline,'/',',');
                        %tline
                        switch i
                            case 1  % line is latitude, longitude in degrees
                                Vlola=sscanf(tline,'%f,%f,%f,%f');
                                LIMITS.lonmin=Vlola(1);  % West
                                LIMITS.lonmax=Vlola(2);  % East
                                LIMITS.latmin=Vlola(3);  % South
                                LIMITS.latmax=Vlola(4);  % North
                            case 2  % line is UTM Easting, Northing in meters
                                Vutm=sscanf(tline,'%f,%f,%f,%f');
                                LIMITS.Emin=Vutm(1);  % West
                                LIMITS.Emax=Vutm(2);  % East
                                LIMITS.Nmin=Vutm(3);  % South
                                LIMITS.Nmax=Vutm(4);  % North
                            case 3  % line is UTM zone
                                LIMITS.UTMzone=str2double(tline);
                            otherwise
                                error('unknown case');
                        end
                    end
                end
            end
        else
            break
        end
    end
    fclose(fid);
else
    error(sprintf('Could not open file named %s\n',fname));
    return
end

