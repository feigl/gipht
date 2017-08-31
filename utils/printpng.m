function printpng(pngfilename)
% print current figure to png file name with 600 DPI

if nargin < 1
    pngfilename = mfilename;
end

if numel(strfind(pngfilename,'.png')) <= 0
    pngfilename = sprintf('%s.png',pngfilename);
end

t1=strrep(pngfilename,'_',' ');
t2=date;
t3=sprintf('%s %s',t1,t2);

% Label with time, date, and file name
% do not do this if suplot has been called
Ch = get(gcf,'Children');
if numel(Ch) == 1
   h=text(1.1,-0.1,t3,...
      'Units','normalized'...
      ,'VerticalAlignment','Bottom'...
      ,'HorizontalAlignment','Left'...
          ,'Interpreter','None'...
      ,'Rotation',90);
end

print (gcf,pngfilename,'-dpng','-r600');

return


