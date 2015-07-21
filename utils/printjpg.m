function printjpg(jpgfilename)
% print current figure to jpg file name with 600 DPI

if nargin < 1
    jpgfilename = mfilename;
end

if numel(strfind(jpgfilename,'.jpg')) <= 0
    jpgfilename = sprintf('%s.jpg',jpgfilename);
end

t1=strrep(jpgfilename,'_',' ');
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


print (gcf,jpgfilename,'-djpeg','-r 600');

return


