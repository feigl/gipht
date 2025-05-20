function plot_and_label_wells(Twells,UTMranges)
% plot Wells on current figure
kilo=1000;

%if isfield(Twells,'Type')
Wtypes={'Production','Injection','Observation','Thermal Gradient'};
Wsymbs={'^','v','o','s'};
Wcolrs={'r','b','g','k'};
for ii=1:numel(Wtypes)
    iType=find(contains(Twells.Type,Wtypes{ii}));
    nTypes= numel(iType)
    if nTypes >  0
        symb1=Wsymbs{ii};
        colr1=Wcolrs{ii};

        % plot all wells
        Hh(ii)=plot(Twells.Easting_m(iType)/kilo,Twells.Northing_m(iType)/kilo ...
            ,'Marker',symb1,'MarkerSize',6,'MarkerFaceColor',colr1...
            ,'LineStyle','none');
        % label wells without letters
        iType=intersect(iType,find(~contains(Twells.Well_ID,'A')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'B')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'C')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'D')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'E')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'F')));
        iType=intersect(iType,find(~contains(Twells.Well_ID,'I')));


        for iii=1:numel(iType)
            for iiii=1:3
                switch iiii
                    case 1
                        xplot=Twells.Easting_m(iType(iii))/kilo;
                        yplot=(UTMranges(3) + 0*iii)/kilo;
                        bg='w';rot=90;
                    case 2
                        xplot=(UTMranges(1) + 0*iii)/kilo;
                        yplot=Twells.Northing_m(iType(iii))/kilo;
                        bg='w';rot=0;
                    case 3
                        xplot=Twells.Easting_m(iType(iii))/kilo;
                        yplot=Twells.Northing_m(iType(iii))/kilo;
                        bg='none';rot=0;
                end
                %Twells.Well_ID(iType(iii))
                text(xplot,yplot,Twells.Well_ID(iType(iii)) ...
                    ,'color','k','BackgroundColor',bg...
                    ,'HorizontalAlignment','center','VerticalAlignment','bottom' ...
                    ,'Interpreter','none','Rotation',rot);
            end
        end
    end
    legend(Hh,Wtypes)
end
% else
%     plot(Twells.Easting_m/kilo,Twells.Northing_m/kilo ...
%                 ,'Marker','*','MarkerSize',6,'MarkerFaceColor','k'...
%                 ,'LineStyle','none');
%     legend('well')
% end
return
end