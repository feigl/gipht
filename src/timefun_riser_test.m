clear all
close all;

% use datetime
% tepochs = datenum(datetime(2017,1,1,12,0,0))+1:100;
% treference1 = datenum(datetime(2017,5,1));
% treference2 = datenum(datetime(2017,6,1));
% use decimal years
tepochs = dyear(2017,1,1)+[1:365]/365;
treference1 = dyear(2017,5,1);
treference2 = dyear(2017,6,1);
ft=timefun_riser(tepochs,treference1,treference2);
figure
plot(tepochs,ft,'or');