function nf = figure4(nf)
%function H = figure4(nf)
% open a square figure window in a quadrant of the full screen
% Based on https://www.mathworks.com/matlabcentral/answers/141266-plotting-two-figures-side-by-side
% 20200504 Kurt Feigl

% get the size of the screen
%screensize = get(0, 'ScreenSize' )
screensize = get( groot, 'ScreenSize' );
%screensize = get(0,'MonitorPositions')
% screensize =
% 
%            1           1        2560        1440
% 
% [left bottom width height]

halfwidth = 0.95*floor(screensize(3)/2);
halfheight = 0.95*floor(screensize(4)/2);
position = screensize;
sidelength = 0.95*min([halfwidth, halfheight]);
position(3) = sidelength;
position(4) = sidelength;

% choose which quadrant
iquad = mod(nf,4);
switch iquad
    case 1 % upper left
        position(1)  = screensize(1)+halfwidth-sidelength;
        position(2)  = screensize(2)+halfheight;
    case 2 % upper right
        position(1) = screensize(1)+halfwidth;
        position(2) = screensize(2)+halfheight;
    case 3 % lower left
        position(1) = screensize(1)+halfwidth-sidelength;
        position(2) = screensize(2);
    case 0 % lower right
        position(1) = screensize(1)+halfwidth;
        position(2) = screensize(2);
    otherwise
        iquad
        warning('invalid iquad');
        position(1) = screensize(1);
        position(2) = screensize(2);
end
%position
H=figure(nf);
%
set(H,'OuterPosition',position);
return
end

