function printpng(pngfilename)
% print current figure to png file name with 600 DPI

if nargin < 1
    pngfilename = mfilename;
end

if numel(strfind(pngfilename,'.png')) <= 0
    pngfilename = sprintf('%s.png',pngfilename);
end

% Label with time, date, and file name
% do not do this if suplot has been called
% Ch = get(gcf,'Children');
% if numel(Ch) == 1
h=text(-0.09,1.09 ...
    ,sprintf('%s %s',pngfilename,datestr(now,30)) ...
    ,'Units','normalized'...
    ,'VerticalAlignment','Top'...
    ,'HorizontalAlignment','Left'...
    ,'Interpreter','None'...
    ,'FontName','Courier','FontSize',9 ...
    ,'Rotation',0);
% end


print (gcf,pngfilename,'-dpng','-r600');

return


