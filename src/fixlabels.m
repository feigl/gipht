function fixlabels(xlab,xfmt,ylab,yfmt,xfontsize,yfontsize,weight)
%function fixlabels(xlabel,xfmt,ylabel,yfmt,xfontsize,yfontsize)
% Print the same number of digits after the decimal point.
%
% x = 1:100;
% y = 5 * (x - mean(x));
% plot(x,y,'ro');
% fixlabels('x label is here ','%.0f','y label is here ','%.1f');
% Kurt Feigl
% 2015-MAR-25 modified to work for MATLAB 8.4.0.150421 (R2014b)

% 2011-NOV-21 Cannot get this to work for AGU
% 20200426 Add optional arguments 5, 6, 7

%       xticklabels=get(gca,'XTickLabel');xtickvals=str2num(xticklabels);
%       for k=1:numel(xtickvals)
%          xticklabels2{k}=sprintf('%9.3f',xtickvals(k));
%       end
%       set(gca,'XTickLabel',xticklabels2);

narginchk(4,7);

if exist('xfontsize','var') == 0
    xfontsize = 10;
end
if exist('yfontsize','var') == 0
    yfontsize = 10;
end
if exist('weight','var') == 0
    weight = 'bold';
end

%if exist('xfmt','var')
if numel(xfmt) > 0
    xlabel(xlab,'FontName','Helvetica','Fontsize',xfontsize,'FontWeight','bold');
    xticklabels=get(gca,'XTickLabel');
%     xtickvals=str2num(xticklabels); % This returns only the mantissa, not the exponent
    xtickvals=get(gca,'XTick'); % 
    for k=1:numel(xtickvals)
        %    xticklabels2{k}=sprintf('%10.4f',xtickvals(k));
        %if mod(k,2) == 0 || numel(xtickvals) < 4
        if mod(k,2) == 0 && k <= numel(xtickvals)
            xticklabels2{k}=sprintf(xfmt,xtickvals(k));
        else
            xticklabels2{k} = sprintf('   ');
        end
        %fprintf(1,'%d %s %f %s\n',k,xticklabels(k),xtickvals(k),char(xticklabels2{k}));
    end
    xticklabels2 = char(xticklabels2);
    %set(gca,'XTickLabelMode','manual');
    
    set(gca,'XTickLabel',xticklabels2);
    set(gca,'FontName','Helvetica','Fontsize',xfontsize,'FontWeight',weight);
else
    xlabel(xlab,'FontName','Helvetica','Fontsize',xfontsize,'FontWeight',weight,'Interpreter','none');
end

%if exist('yfmt','var')
if numel(yfmt) > 0
    ylabel(ylab,'FontName','Helvetica','Fontsize',yfontsize,'FontWeight',weight);
    yticklabels = get(gca,'YTickLabel');
    %ytickvals=str2num(yticklabels) % This returns only the mantissa, not the exponent
    ytickvals = get(gca,'YTick');
    for k=1:numel(ytickvals)
        %yticklabels2{k}=sprintf('%10.4f',ytickvals(k));
        %ytickvals(k)=char(str2num(yticklabels{k}));
        %if mod(k,2) == 1  || numel(ytickvals) < 4
         if mod(k,2) == 1  && k <= numel(ytickvals);
            ytlab1 = sprintf(yfmt,ytickvals(k));
            ytlab1 = strrep(ytlab1,'E','\times 10^{');
            if numel(strfind(ytlab1,'{')) > 0
                ytlab1 = strcat(ytlab1,'}');
            end
            yticklabels2{k}=ytlab1;
        else
            yticklabels2{k} = sprintf('    ');
        end
%         fprintf(1,'%d %s %f %s\n',k,char(yticklabels{k}),ytickvals(k),char(yticklabels2{k}));
    end
    %set(gca,'YTickLabelMode','manual');
    yticklabels2 = char(yticklabels2);
    set(gca,'YTickLabel',yticklabels2);
    set(gca,'FontName','Helvetica','Fontsize',yfontsize,'FontWeight',weight);

    %h = gca;
% else
%     set(gca,'YTickLabelMode','auto');
else
    ylabel(ylab,'FontName','Helvetica','Fontsize',yfontsize,'FontWeight',weight,'Interpreter','none');
end

%h = gca;
%h = 1;

return

end

