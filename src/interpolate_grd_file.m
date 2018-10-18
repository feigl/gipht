function output_grd_file_name = interpolate_grd_file(input_grd_file_name,xax2,yax2)
%% calculate UTM coordinates and write a new grid file
% both input and output are GMT grid files
% 20180926 Kurt Feigl

idebug = 0;


%% check 
INFO = grdinfo3(input_grd_file_name);


%% read file 
[xax1,yax1,zvals1]=grdread3(input_grd_file_name);
[GX1,GY1]=meshgrid(xax1,yax1);
[nrows1,mcols1] = size(GX1);

%% make new mesh
% nrows2 = numel(xax2);
% mcols2 = numel(yax2);
[GX2,GY2]=meshgrid(xax2,yax2);
[nrows2,mcols2] = size(GX2);

%% handle phase differently, as real and imaginary parts
if contains(INFO.zname,'phase','IgnoreCase',true) == 1  && contains(INFO.zname,'unwrapped','IgnoreCase',true) == 0
    %whos
    % construct interpolant functions
    Fx=scatteredInterpolant(colvec(GX1),colvec(GY1),colvec(sin(zvals1)),'linear','none');   
    Fy=scatteredInterpolant(colvec(GX1),colvec(GY1),colvec(cos(zvals1)),'linear','none');
    % perform the interpolation
    xvals = Fx(colvec(GX2),colvec(GY2));
    yvals = Fy(colvec(GX2),colvec(GY2));
    % reshape
    zvals2 = reshape(angle(complex(xvals,yvals)),nrows2,mcols2);
else
    Fz=scatteredInterpolant(colvec(GX1),colvec(GY1),colvec(zvals1),'linear','none');
    zvals2 = Fz(colvec(GX2),colvec(GY2));
    zvals2 = reshape(zvals2,nrows2,mcols2);
end

%% write new header
INFO.dx =  xax2(2)-xax2(1)
INFO.dy =  yax2(2)-yax2(1)
INFO.xmin =  nanmin(xax2);
INFO.xmax =  nanmax(xax2);
INFO.ymin =  nanmin(yax2);
INFO.ymax =  nanmax(yax2);
INFO.zmin =  nanmin(colvec(zvals2));
INFO.zmax =  nanmax(colvec(zvals2));
INFO.nx =    mcols2;
INFO.ny =    nrows2;
INFO.xname = 'Easting in meters';
INFO.yname = 'Northing in meters';
INFO.ispixelreg = 1;

if idebug == 1
    figure;hold on;
    imagesc(xax1,yax1,zvals1);
    axis xy
    axis image
    axis equal
    
    title('before interpolation');
    colormap(cmapgraynan);
    colorbar;
    
    
    figure;
    imagesc(xax2,yax2,zvals2);
    axis xy
    axis image
    axis equal
    colormap(cmapgraynan);
    colorbar;
    title('after interpolation');
end



%% Write a GMT grid file
output_grd_file_name = strrep(input_grd_file_name,'.grd',sprintf('_interp%05dx%05d.grd',nrows2,mcols2));
grdwrite3(xax2,yax2,zvals2,output_grd_file_name,INFO);
fprintf(1,'Wrote %s\n',output_grd_file_name);
dirlist = dir(output_grd_file_name);

%% Make a map of it
if idebug == 1
   H=map_grd(output_grd_file_name,cmapgraynan);
end


return
end

