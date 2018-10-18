function output_grd_file_name = cut_grd_file(input_grd_file_name,xmin,xmax,ymin,ymax)
%% cut a grid and write a new grid file
% both input and output are GMT grid files
% 20180926 Kurt Feigl

%% 
INFO = grdinfo3(input_grd_file_name)

%% read file 
[xax1,yax1,zvals1]=grdread3(input_grd_file_name);
[nrows1,mcols1] = size(zvals1)
[XAX1,YAX1] = meshgrid(xax1,yax1);



%% make new axes
xmin2=xmin-mod(xmin,INFO.dx);
xmax2=xmax+mod(xmax,INFO.dx);
ymin2=ymin-mod(ymin,INFO.dy);
ymax2=ymax+mod(ymax,INFO.dy);
xax2 = xmin2:INFO.dx:xmax2;
yax2 = ymin2:INFO.dy:ymax2;

%% cut
ii = intersect(find(xax1 >= xmin2),find(xax1 <= xmax2));
jj = intersect(find(yax1 >= ymin2),find(yax1 <= ymax2));
zvals2 = zvals1(jj,ii);
[nrows2,mcols2] = size(zvals2)

% figure;hold on;
% imagesc(xax1,yax1,zvals1);
% axis xy
% axis image
% axis equal
% 
% title('uncut');
% colormap(jet);
% colorbar;
% 
% 
% figure;
% imagesc(xax2,yax2,zvals2);
% axis xy
% axis image
% axis equal
% colormap(jet);
% colorbar;
% title('cut');


%% write new header
INFO.xmin =  xmin2;
INFO.xmax =  xmax2;
INFO.ymin =  ymin2;
INFO.ymax =  ymax2;
INFO.zmin =  nanmin(colvec(zvals2));
INFO.zmax =  nanmax(colvec(zvals2));
INFO.nx =    mcols2;
INFO.ny =    nrows2;


%% Write a GMT grid file
output_grd_file_name = strrep(input_grd_file_name,'.grd',sprintf('_cut__%05dx%05d.grd',nrows2,mcols2));
grdwrite3(xax2,yax2,zvals2,output_grd_file_name,INFO);
fprintf(1,'Wrote %s\n',output_grd_file_name);
dirlist = dir(output_grd_file_name)

%% Make a map of it
%H=map_grd(output_grd_file_name,jet);


    

return
end

