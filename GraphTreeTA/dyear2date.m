function [ Yr, M, D, DOY] = dyear2date(Dy)
% function [ Yr, M, D, DOY] = dyear2date(Dy)
% Converts decimal years found with dyear.m back into calendar dates; given
% in vector formats
% 
% Elena C. Baluyut,  UW-Madison
% 2015-02-09

% Find number of dates to convert
x = numel(Dy);

% Loop through each date
for i = 1:x
    dy = Dy(i);
    yr = floor(dy); %finds current calendar year
    yr_next = ceil(dy); %finds next calendar year
    
    frac = dy - yr; %finds fractional portion of year
    doy = ((frac.*(datenum([yr_next 0 0 0 0 0])-datenum([yr 0 0 0 0 0])))+1); %finds day of year in accordance with date2doy script line 45
    
    
    if mod(yr, 4) == 0  %determine if year is leap year
        lp_yr = 1;
    else
        lp_yr = 0;
    end
    
    d_markers = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]; %day count for end of each month
    
    if isequal(lp_yr, 1) == 1
        d_markers(2:end) = d_markers(2:end) + lp_yr; %adjust for extra day in February
    end
    
    dom = find(d_markers <= doy, 1, 'last'); %find day of year at beginning of the month
    
    if isempty(dom) == 1
        m = 1;
        d = round(doy); %for days in January
    else
        m = dom+1;
        d = round(doy-d_markers(dom)); %for all others, find day of month
    end
    
    % store data
    Yr(i,1) = yr;
    M(i,1) = m;
    D(i,1) = d;
    DOY(i,1) = doy;
end

return

