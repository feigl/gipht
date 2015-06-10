function h = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab)
%function h = draw_profile(xt,yt1,yt2,xlab,y1lab,y2lab,tlab)
% Draw a profile


if numel(xt) < 2 || numel(yt1) < 2 || numel(yt2) < 2
    error('Need at least 2 points to plot a profile');
end


% absolute value of residual
yt3 = abs(yt2-yt1);

h=figure;
subplot(3,1,1);
title(tlab,'FontName','Times','Fontsize',9);
axis off

subplot(3,1,2);
iok = find(isfinite(yt1)==1);
iok = intersect(iok,find(isfinite(yt2)==1));
iok = intersect(iok,find(abs(yt1)>0.0));
iok = intersect(iok,find(abs(yt2)>0.0));

if numel(iok) < 1
    warning('too few points','numel(iok) = %d LT 1\n',numel(iok));
    return;
end
xt=xt(iok);
yt1=yt1(iok);
yt2=yt2(iok);
yt3=yt3(iok);

plot(xt,yt1,'ro','MarkerFaceColor','r','MarkerSize',4);
axis([nanmin(nanmin(xt)) nanmax(nanmax(xt)) -Inf +Inf]);
hold on;
if numel(strfind(lower(y1lab),'range')) > 0
    axis ij; % make uplift upwards
else
    axis xy; % make displacement positive upwards
end
plot(xt,yt2,'k-','LineWidth',2);
set(gca,'FontName','Helvetica-Bold','Fontsize',12);
%legend('Observed','Modeled','Location','NorthEast');
legend('Observed','Modeled','Location','Best');
xlabel(xlab,'FontName','Helvetica-Bold','Fontsize',12);
%ylabel(y1lab,'FontName','Helvetica-Bold','Fontsize',12);
title(y1lab,'FontName','Helvetica-Bold','Fontsize',12);
ylabel(y2lab,'FontName','Helvetica-Bold','Fontsize',12);
%    fixlabels('Easting (km)','%4.0f','Range change (mm)','%.0f');

subplot(3,1,3);
plot(xt,yt3,'b.','LineWidth',2);
axis([min(xt) max(xt) -Inf +Inf]);hold on;

%axis ij; % make uplift upwards
set(gca,'FontName','Helvetica-Bold','Fontsize',12,'FontWeight','bold');
xlabel(xlab,'FontName','Helvetica-Bold','Fontsize',12)
ylabel(y2lab,'FontName','Helvetica-Bold','Fontsize',12)
%fixlabels('Easting (km)','%4.0f','Range change (mm)','');
%   set(gca,'FontName','Helvetica-Bold','Fontsize',12,'FontWeight','bold');
%legend('Angular Deviation','Wrapped Residual','Location','NorthEast');
%legend('Absolute Deviation','Location','NorthEast');
legend('Absolute Deviation','Location','Best');

return

end

