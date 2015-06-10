function fixlabels(xlab,xfmt,ylab,yfmt)
%function h = fixlabels(xlabel,xfmt,ylabel,yfmt)
% Print the same number of digits after the decimal point.
%
% x = 1:100;
% y = 5 * (x - mean(x));
% plot(x,y,'ro');
% fixlabels('x label is here ','%.0f','y label is here ','%.1f');
% Kurt Feigl
% 2015-MAR-25 modified to work for MATLAB 8.4.0.150421 (R2014b)

% 2011-NOV-21 Cannot get this to work for AGU

%       xticklabels=get(gca,'XTickLabel');xtickvals=str2num(xticklabels);
%       for k=1:numel(xtickvals)
%          xticklabels2{k}=sprintf('%9.3f',xtickvals(k));
%       end
%       set(gca,'XTickLabel',xticklabels2);

%if exist('xfmt','var')
if numel(xfmt) > 0
    xlabel(xlab,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
    xticklabels=get(gca,'XTickLabel');
%     xtickvals=str2num(xticklabels); % This returns only the mantissa, not the exponent
    xtickvals=get(gca,'XTick'); % 
    for k=1:numel(xtickvals)
        %    xticklabels2{k}=sprintf('%10.4f',xtickvals(k));
        if mod(k,2) == 0 || numel(xtickvals) < 4
            xticklabels2{k}=sprintf(xfmt,xtickvals(k));
        else
            xticklabels2{k} = sprintf('   ');
        end
        %fprintf(1,'%d %s %f %s\n',k,xticklabels(k),xtickvals(k),char(xticklabels2{k}));
    end
    xticklabels2 = char(xticklabels2);
    %set(gca,'XTickLabelMode','manual');
    
    set(gca,'XTickLabel',xticklabels2);
    set(gca,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
else
    xlabel(xlab,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
end

%if exist('yfmt','var')
if numel(yfmt) > 0
    ylabel(ylab,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
    yticklabels=get(gca,'YTickLabel');
    %ytickvals=str2num(yticklabels) % This returns only the mantissa, not the exponent
    ytickvals  =get(gca,'YTick');
    for k=1:numel(yticklabels)
        %yticklabels2{k}=sprintf('%10.4f',ytickvals(k));
        %ytickvals(k)=char(str2num(yticklabels{k}));
         if mod(k,2) == 1  || numel(ytickvals) < 4
            ytlab1 = sprintf(yfmt,ytickvals(k));
            ytlab1 = strrep(ytlab1,'E','\times 10^{');
            if numel(strfind(ytlab1,'{')) > 0
                ytlab1 = strcat(ytlab1,'}');
            end
            yticklabels2{k}=ytlab1;
        else
            yticklabels2{k} = sprintf('    ');
        end
        fprintf(1,'%d %s %f %s\n',k,char(yticklabels{k}),ytickvals(k),char(yticklabels2{k}));
    end
    %set(gca,'YTickLabelMode','manual');
    yticklabels2 = char(yticklabels2);
    set(gca,'YTickLabel',yticklabels2);
    set(gca,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
    %h = gca;
% else
%     set(gca,'YTickLabelMode','auto');
else
    ylabel(ylab,'FontName','Helvetica','Fontsize',10,'FontWeight','bold');
end

%h = gca;
%h = 1;

return

end

