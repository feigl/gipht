function printpdf(pdffilename,orientation)
%function printpdf(pdffilename,orientation)
% write current graphics window to a PDF file
% last update: 20130520 Kurt

if nargin < 1
    pdffilename = mfilename;
end

if nargin ==2 && exist('orientation','var') == 1
        if strcmp(lower(orientation),'landscape')== 1
        %% force landscape
        set(gcf,'PaperUnits','inch');
        set(gcf,'PaperPosition', [0. 0. 11 8.5]);
        set(gcf,'PaperOrientation', 'landscape');
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperSize',[11.0 8.5]);
    end
end


if numel(strfind(pdffilename,'.pdf')) <= 0
    pdffilename = sprintf('%s.pdf',pdffilename);
    % 20130514 Kurt: Change slashes to back slashes if needed
    pdffilename = strrep(pdffilename,'/',filesep);
    % remove blanks from file name
    pdffilename = strrep(pdffilename,' ','');
end

% print current figure to pdf file name with 1200 DPI and TIFF
t0=pwd;
%  20130514 Kurt: this is the title. Replace the underscore to avoid Tex errors
%t1=strrep(pdffilename,'_','\_');
% 20130620 - replacement is done again below. Once is enough
t1=pdffilename;
%t2=date;
t2=datestr(now,31); %31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
tu=getenv('USER');
%t3=sprintf('%s %s %s %s',t1,t0,t2,tu);
%t3=sprintf('%s %s %s\n%s',t1,t2,tu,t0);
t3=sprintf('%s %s %s',t1,t2,tu);
%t4= strrep(t3,'_',[filesep '_']);
% 20130514 Kurt: Replace the underscore to avoid Tex errors
%t4= strrep(t3,'_','\_');


%% label with file name, date and time on bottom
subplot('Position',[0 0 20 0.05],'Units','Centimeters');
axis off
text(10,0,t3...
    ,'Units','Centimeters'...
    ,'VerticalAlignment','Bottom'...
    ,'HorizontalAlignment','Center'...
    ,'Clipping','off'...
    ,'FontName','Courier','FontSize',9 ...
    ,'Rotation',0 ...
    ,'Interpreter','None');
% %% label with file name, date and time on top
% subplot('Position',[0.01 0.95 0.99 0.05],'Units','Normalized');
% axis off
% text(0,0.,t3...
%     ,'Units','Normalized'...
%     ,'VerticalAlignment','Bottom'...
%     ,'HorizontalAlignment','Left'...
%     ,'Clipping','off'...
%     ,'FontName','Courier','FontSize',9 ...
%     ,'Rotation',0 ...
%     ,'Interpreter','None');

% %% label with file name, date and time on side
% subplot('Position',[0.01 0.01 0.05 0.99],'Units','Normalized');
% axis off
% % coordinates for text are inside the rectangle defined by subplot above
% text(0.1,0.5,t3...
%     ,'Units','Normalized'...
%     ,'VerticalAlignment','Top'...
%     ,'HorizontalAlignment','Center'...
%     ,'Clipping','off'...
%     ,'FontName','Courier','FontSize',9 ...
%     ,'Rotation',90 ...
%     ,'Interpreter','None');

% Label the figure
% Does not work on Hengill
% mycomputer = computer;
% if strcmp(mycomputer, 'GLNXA64')==0
% http://www.mathworks.com/support/solutions/en/data/1-703XB9/?solution=1-703XB9
% Subject:
% 
% Is it possible to programmatically check whether MATLAB has been started
% with the "-nodisplay" option? Problem Description:
% 
% I have a program that needs to behave differently depending on whether
% MATLAB has a display or not. However, I cannot determine a way to
% programmatically check this.
% 
% The reason I need to do this is that I sometimes start MATLAB in batch
% mode from a shell script, for testing my programs.
% 
% Solution:
% 
% The ability to programmatically check whether MATLAB has a display is not
% available in MATLAB 7.6 (R2008a).
% 
% As a workaround, you can do one of the following:
% 
% 1. Manually set or unset an environment variable in your shell script
% that launches MATLAB, so you can query it from inside MATLAB. For
% instance, setting the ISDISPLAY environment variable in your shell script
% before launching MATLAB:
% 
% 
% setenv ISDISPLAY no matlab -r "foo; quit"
% 
% (here shown with C Shell syntax), means that the MATLAB command
% getenv('ISDISPLAY') will return the string 'no' in any MATLAB processes
% that are spawned from this shell.
% 
% 2. Query the "ScreenSize" property of the root object inside MATLAB:
% get(0, 'ScreenSize')
% When there is no display, this returns [1 1 1 1] instead of an actual screen size. However, this relies on behavior that isn't actually specified (by the doc, for instance) to work in any particular way, so may be subject to change in the future. If you were going to use this many times, it might be wise to wrap it in a function (e.g. create an "isdisplay.m" function file), so you can easily change the implementation in the future, if needed. (This method worked as of MATLAB R2008a.) 
% ss4 = get(0, 'ScreenSize');
% 
% 
% if ss4 == [1 1 1 1]
%    if ismac == 1
%        figfilename = strrep(pdffilename,'pdf','fig');
%        saveas(gcf,figfilename,'fig');
%    else
%    psfilename = strrep(pdffilename,'pdf','ps');
%    fprintf(1,'Printing PostScript to file named %s\n',psfilename);
%    print(psfilename,'-dpsc2','-r1200'); % print PS if no display  
% else
%   print(pdffilename,'-dpdf','-r1200')
% end


try
    fprintf(1,'Printing PDF to file named %s\n',pdffilename);   
    print(gcf,pdffilename,'-dpdf','-r1200'); % print PDF   
catch
    warning(sprintf('ERROR in %s. Trying to catch...\n',mfilename));   
    figfilename = strrep(pdffilename,'pdf','fig');
    saveas(gcf,figfilename,'fig');
end


%     try
%         psfilename = strrep(pdffilename,'pdf','ps');
%         fprintf(1,'Printing PostScript to file named %s\n',psfilename);
%         print(psfilename,'-dpsc2','-r1200'); % print PS if no display
%     catch
%         warning(sprintf('ERROR in %s\n',mfilename));
%         return
% %         jpgfilename = strrep(psfilename,'ps','jpg');
% %         fprintf(1,'Printing JPEG to file named %s\n',jpgfilename);
% %         print(jpgfilename,'-djpeg','-r1200'); % print JPG    
%     end

   
% end

% if exist('ghandle','var') == 1
%     print(ghandle,pdffilename,'formattype','dpdf','resolution',1200); % otherwise, print PDF
% else
%     %print(gcf,pdffilename,'-dpdf','-r1200'); % otherwise, print PDF
%     %print(pdffilename,'-dpdf','-r1200'); % otherwise, print PDF
%     print(gcf,pdffilename,'formattype','dpdf','resolution',1200); % otherwise, print PDF   
% end

%print(pdffilename,'formattype','dpdf','resolution',1200); % otherwise, print PDF   
%print(gcf,pdffilename,'-dpdf','-r1200'); % otherwise, print PDF   

return


