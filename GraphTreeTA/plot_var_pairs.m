function [ minWd, h] = plot_var_pairs(Wd, tm, ts, ylabeltxt, titlestr)
% function [ minWd, hbar] = plot_var_pairs(Q, V, tu)
% Plots the relative variance for each epoch as a bar graph with epochs labeled with
% corresponding calendar dates
% 
% INPUTS:
%   Wd - vector of variances for pair-wise measurements (to plot st. dev.,
%        enter sqrt(Wd))
%   tu - vector of epochs
%   ylabeltxt - string containing label for y axis
%   titlestr - string containing title
%
% OUTPUTS:
%   minWd - minimum uncertainty of pair-wise measurements (used for
%           normalizing)
%   hbar  - handle for bar graph 
%
% Elena C. Baluyut, UW-Madison
% 2015-07-21

% Calculate variance of epoch-wise measurements


minWd = min(Wd); % locate smallest value to normalize by
Wd_norm = Wd./minWd; % normalize standard deviations

index = 1:numel(tm); % assign epoch index numbers

% Display 
figure;
h = bar(index, Wd_norm);
%set(h, 'Interpreter', 'tex');
%ht = title('Standard Deviation in Epoch-wise Measurements by Epoch Index')
ht = title(titlestr);
hx = xlabel('pair index');
%hy = ylabel(sprintf('uncertainty in %s', yunit));
hy = ylabel(ylabeltxt);
set(hy, 'Interpreter', 'tex');
set(ht,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
set(hy,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
set(hx,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');

% Find calendar dates from decimal years 

    [ cal_string_tm ] = print_caldate( tm );
    [ cal_string_ts ] = print_caldate( ts );
    cal_string = strcat(cal_string_tm, ' to ', cal_string_ts);

% Print calendar dates to figure above each bar
for i = 1:numel(tm)
    hText = text(index(i)+.55, Wd_norm(i),  sprintf('  %s', cal_string{i}));
    set(hText,'VerticalAlignment','bottom', 'rotation', 90, 'FontSize',10); 
end

% Extend figure to fit text 
text_pix = 2*ceil(4/3*10); % calculate pixel size of string 
[ax] = axis;
axpx_ratio = ax(4)/80; % ratio of y axis units to pixels
maxbar_pix = 80*max(Wd_norm)/ax(4); % find pixel length of longest bar
axis_pix = maxbar_pix+text_pix; % pixel length
axis_len = axis_pix*axpx_ratio; % convert to axis length
axis([min(index)-.5, max(index)+1, 0, axis_len]);

% return graphics handle to current figure;
h = gcf;

end

