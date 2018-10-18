function output_grd_file_name = project_grd_to_utm(input_grd_file_name)
%% calculate UTM coordinates and write a new grid file
% both input and output are GMT grid files
% 20180926 Kurt Feigl

%% check that coordinates are geographic
INFO = grdinfo3(input_grd_file_name);

if contains(INFO.xname,'Longitude','IgnoreCase',true) ~= 1 || contains(INFO.yname,'Latitude','IgnoreCase',true) ~= 1
    error(('Input grid files should be in geographic coordinates (latitude, longitude) in degrees'));
end

%% read file and make grid in geographic coordinates
[lonax,latax,zvals]=grdread3(input_grd_file_name);
[ngrows,mgcols] = size(zvals);
[LON,LAT] = meshgrid(lonax,latax);


%% make regular grid in UTM coordinates
[Eutm,Nutm,UTMzone] = deg2utm(colvec(LAT),colvec(LON));
dE = (nanmax(colvec(Eutm))-nanmin(colvec(Eutm)))/mgcols
dN = (nanmax(colvec(Nutm))-nanmin(colvec(Nutm)))/ngrows

E1 = nanmin(colvec(Eutm));
E2 = nanmax(colvec(Eutm));
N1 = nanmin(colvec(Nutm));
N2 = nanmax(colvec(Nutm));
Eax = E1:dE:E2;
Nax = N1:dN:N2;
murows = numel(Eax)
nucols = numel(Nax)


[GE,GN]=meshgrid(Eax,Nax);
[nurows,mucols] = size(GE);


%% handle phase differently, as real and imaginary parts
if contains(INFO.zname,'phase','IgnoreCase',true) == 1  && contains(INFO.zname,'unwrapped','IgnoreCase',true) == 0
    whos
    % construct interpolant functions
    Fx=scatteredInterpolant(colvec(Eutm),colvec(Nutm),colvec(sin(zvals)),'linear','none');   
    Fy=scatteredInterpolant(colvec(Eutm),colvec(Nutm),colvec(cos(zvals)),'linear','none');
    % perform the interpolation
    xvals = Fx(colvec(GE),colvec(GN));
    yvals = Fy(colvec(GE),colvec(GN));
    % reshape
    zvals = reshape(angle(complex(xvals,yvals)),nurows,mucols);
else
    Fz=scatteredInterpolant(colvec(Eutm),colvec(Nutm),colvec(zvals),'linear','none');
    zvals = Fz(colvec(GE),colvec(GN));
    zvals = reshape(zvals,nurows,mucols);
end

%% write new header
INFO.dx =  dE;
INFO.dy =  dN;
INFO.xmin =  nanmin(Eax);
INFO.xmax =  nanmax(Eax);
INFO.ymin =  nanmin(Nax);
INFO.ymax =  nanmax(Nax);
INFO.zmin =  nanmin(colvec(zvals));
INFO.zmax =  nanmax(colvec(zvals));
INFO.nx =    mucols;
INFO.ny =    nurows;
INFO.xname = 'Easting in meters';
INFO.yname = 'Northing in meters';


%% Write a GMT grid file
output_grd_file_name = strrep(input_grd_file_name,'.grd','_utm.grd');
grdwrite3(Eax,Nax,zvals,output_grd_file_name,INFO);
fprintf(1,'Wrote %s\n',output_grd_file_name);
dirlist = dir(output_grd_file_name)

%% Make a map of it
%H=map_grd(output_grd_file_name,cmapgraynan);


    

return
end

