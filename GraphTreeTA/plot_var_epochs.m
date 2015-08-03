function [ minWd, h] = plot_var_epochs(Wd, tu)
% function [ minWd, hbar] = var_epochs(Q, V, tu)
% Plots the relative variance for each epoch as a bar graph with epochs labeled with
% corresponding calendar dates
% 
% INPUTS:
%   Wd - vector of variances for epoch-wise measurements (to plot st. dev.,
%        enter sqrt(Wd))
%   tu - vector of epochs
%
% OUTPUTS:
%   minWd - minimum uncertainty of epoch-wise measurements (used for
%           normalizing)
%   hbar  - handle for bar graph 
%
% Elena C. Baluyut, UW-Madison
% 2015-07-21

% Calculate variance of epoch-wise measurements


minWd = min(Wd); % locate smallest value to normalize by
Wd_norm = Wd./minWd; % normalize standard deviations

index = 1:numel(tu); % assign epoch index numbers

% Display 
figure;
h = bar(index, Wd_norm);
title('Variance in Epoch-wise Measurements by Epoch Index')
xlabel('epoch index')
h = ylabel('uncertainty in \DeltaV');
set(h, 'Interpreter', 'tex');

% Find calendar dates from decimal years
[ cal_string ] = print_caldate( tu );

% Print calendar dates to figure above each bar
for i = 1:numel(tu)
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

