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
% mcols2 = numel(xax2);
% nrows2 = numel(yax2);
% [XAX2,YAX2] = meshgrid(xax2,yax2);

%% cut
% II = intersect(find(XAX1 >= xmin2),find(XAX1 <= xmax2));
% JJ = intersect(find(YAX1 >= ymin2),find(YAX1 <= ymax2));
% whos 
% zvals2 = zvals1(JJ,II);
%zvals2 = zvals1(intersect(II,JJ));

ii = intersect(find(xax1 >= xmin2),find(xax1 <= xmax2));
jj = intersect(find(yax1 >= ymin2),find(yax1 <= ymax2));
zvals2 = zvals1(jj,ii);

%[nrows2,mcols2] = size(zvals2)

% ii1 = floor((xmin2-min(xax1))/INFO.dx)
% ii2 =  ceil((max(xax1)-xmax2)/INFO.dx)
% jj1 = floor((ymin2-min(yax1))/INFO.dy)
% jj2 =  ceil((max(yax1)-ymax2)/INFO.dy)
% zvals2 = zvals1(jj1:jj2,ii1:ii2);
[nrows2,mcols2] = size(zvals2)


figure;hold on;
imagesc(xax1,yax1,zvals1);
axis xy
axis image
axis equal

%plot(meshgrid(xax1,yax1),'r*-');
plot(XAX2,YAX2,'k.-');
title('uncut');
colormap(jet);
colorbar;


figure;
imagesc(xax2,yax2,zvals2);
axis xy
axis image
axis equal
colormap(jet);
colorbar;
title('cut');


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
output_grd_file_name = strrep(input_grd_file_name,'.grd',sprintf('_%05dx%05d.grd',nrows2,mcols2));
grdwrite3(xax2,yax2,zvals2,output_grd_file_name,INFO);
fprintf(1,'Wrote %s\n',output_grd_file_name);
dirlist = dir(output_grd_file_name)

%% Make a map of it
H=map_grd(output_grd_file_name,jet);


    

return
end

