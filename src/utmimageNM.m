function utmimageNM(ims,tls,nrows,mcols....
            ,wesn,titlestr,climit,dotutmx,dotutmy,ctab)
% function utmmageNM(ims,tls,nrows,mcols....
%             ,wesn,titlestr,climit,dotutmx,dotutmy,ctab)
%
%          make multi-panel set of UTM maps
%
% Example
%          im=randn([200, 300]);
%          for j=1:6;ims(j,:,:)=im;end
%          size(ims)
%          utmimageNM(ims,{'A','B','C','E','E','F'},2,3,[30.1 30.2 40.1 40.15]...
%              ,'title goes here',[-Inf Inf],0,0,'jet')
%
% Kurt Feigl 
% 2008 MAY 05 original
% 2011 JUL 03 plot top to bottom
% 2012 JUL 05 adapt for rectangular plots

[nimg,ndum,mdum] = size(ims);

if nimg ~= nrows * mcols
   error 'Misount number of images'
end

xutmmin = wesn(1);
xutmmax = wesn(2);
yutmmin = wesn(3);
yutmmax = wesn(4);

yutmdif = yutmmax-yutmmin;
xutmdif = xutmmax-xutmmin;
rat = (8/11)*(yutmdif/xutmdif);
%rat = 1;

nodot = 100;

h1=figure;hold on;
%set(h1,'PaperPosition', [0.25 0.25 8 8]);
set(h1,'PaperPosition', [0.25 0.25 8 11]);

set(h1,'PaperPositionMode','manual');
%set(h1,'PaperSize',[8.5 8.5]);
set(h1,'PaperSize',[8.5 11.0]);


set(h1,'Position',[1 1 1000 1000*rat]);
%set(h1,'Position',[1 1 1000 1000]);


set(h1,'PaperUnits','centimeters');
set(h1,'Units','centimeters');
% set(h1,'PaperUnits','inches');
% set(h1,'Units','inches');


%set(h1,'PaperPosition', [0.25 0.25 8 8*rat]);
wide = 20; %cm
%high = 27; %cm 
high = wide * rat;
while high > 27 || wide > 20
   wide = 0.9 * wide;
   high = wide * rat;
end
% wide = 8; 
% high = 10; % inches 
set(h1,'PaperPosition', [0.5 0.5 wide high]); % LL and UR corners in cm

% noe work in normalized units
wide = 1;
high = 1; 
rat = 1;


set(h1,'PaperPositionMode','manual');
set(h1,'PaperType','usletter');
%set(h1,'PaperSize',[8.5 8.5]);

%set(h1,'Position',[1 1 1000 1000*rat]);
%set(h1,'Position',[1 1 1000 1000*rat]);

% dwide = 1.0/(nrows+1);
% dhigh = rat/(mcols+1);


k=0;
im1 = zeros(nrows*ndum,mcols*mdum);

for i=1:nrows % plot from top to bottom
    for j=1:mcols % plot from left to right
        k=k+1;
        im1((i-1)*ndum+1:i*ndum,(j-1)*mdum+1:j*mdum)=ims(k,:,:);
    end
end
imagesc(im1);axis ij; axis tight; axis equal;
title(titlestr);
xlabel('pixel index for columns');
ylabel('pixel index for rows');
colorbar;

return;





