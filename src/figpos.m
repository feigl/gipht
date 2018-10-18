function H = figpos(ii,nn)
% return position for figure i of n
%position0 = get(0,'MonitorPositions')
pos0 = get(0,'ScreenSize');
nrows = 2;
mcols = ceil(nn/2);

iii=ii;
if iii > nn
    iii = mod(iii,nn);
end
if iii == 0
    iii = 1;
end

cc = 0.99;
width  = floor(cc * pos0(3)/mcols);
height = floor(cc * pos0(4)/nrows);

%[irow,jcol] = ind2sub([nrows,mcols],mod(i,n)+1);

% irow = mod(i,mcols);
% jcol = mod(i,nrows);
A=reshape(1:nn,[mcols,nrows])';
[irow,jcol]=find(A==iii);
irow = irow(1);
jcol = jcol(1);
left   = pos0(1) +          floor((jcol-1)*width/cc);
bottom = pos0(2) + pos0(4) - floor((irow)*height/cc);
H=figure(ii);
pos1 = [left bottom width height]
if numel(pos1) == 4
    fprintf(1,'%5d %5d %5d %5d %5d %5d %5d\n',ii,irow,jcol,pos1);
    set(H,'OuterPosition',pos1,'ToolBar','none','MenuBar','none');
    %,'NumberTitle','off'
    % 'Outerposition
end
return
end

