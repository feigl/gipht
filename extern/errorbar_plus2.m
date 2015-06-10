function h = errobar_plus2(x,y,dx,dy,prop,MarkerSize,LineWidth)
%function h = errobar_plus2(x,y,dx,dy,prop,MarkerSize,LineWidth)
% By Nguyen Chuong 2004/06/30
% modified 2011-06-13 to return a graphics handle - Kurt Feigl
% This is to the problem of legend in errorbar
% This function puts an errobar range on plot

nargchk(3,6,nargin);

% make all data column vectors
x = reshape(x,numel(x),1);
y = reshape(y,numel(y),1);
dx = reshape(dx,numel(dx),1);
dy = reshape(dy,numel(dy),1);

if size(x) ~= size(y)
    error('Dimensions are not consistent!');
end
if ~exist('dy','var')
    dy = dx;
    dx = zeros(size(y));
end
if numel(dx) == 1
    dx = dx * ones(size(x));
end
if numel(dy) == 1
    dy = dy * ones(size(y));
end
if size(dx) ~= size(x)
    error('Dimensions are not consistent!');
end    
if size(dy) ~= size(y)
    error('Dimensions are not consistent!');
end    
if ~exist('MarkerSize','var')
    MarkerSize = 3;
end
if ~exist('LineWidth','var')
    LineWidth = 1;
end
if exist('prop')
	color_list='bgrcmyk';
	mark_list='.ox+*sdv^<>ph';
	color='b';
    for j=1:numel(prop)
    	for i=1:numel(color_list)
            if prop(j)==color_list(i)
               color=color_list(i);
               break
            end
        end
        if prop(j)==color_list(i)
           break
        end
	end
	mark='o';
    for j=1:numel(prop)
    	for i=1:numel(mark_list)
            if prop(j)==mark_list(i)
               mark=mark_list(i);
               break
            end
        end
        if prop(j)==mark_list(i)
           break
        end
	end
else
	color='b'; %default values
	mark='o';
end

hold on


%for i=1:min([numel(x) numel(y) numel(dx) numel(dy)])
for i=1:numel(x)
    h = plot([x(i)-dx(i) x(i)+dx(i)], [y(i)       y(i)      ],color);
    set(h,'LineWidth',LineWidth);
    h = plot([x(i)       x(i)      ], [y(i)-dy(i) y(i)+dy(i)],color);
    set(h,'LineWidth',LineWidth);
end

% draw symbol at midpoint
h = plot(x,y,prop);
set(h,'MarkerSize',MarkerSize);
set(h,'MarkerFaceColor',color);    

%%%% DO NOT DO THIS!
%%%hold off
return
