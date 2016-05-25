function monthnum = monthstr2monthnum(monthstr)
%function monthnum = monthstr2monthnum(monthstr)
% convert month in 3-letter text string in English or French to integer number
% 20160525 Kurt Feigl

% remove insignificant white space
monthstr = strtrim(monthstr);

% truncate to 4 letters
if numel(monthstr) > 3
    monthstr = monthstr(1:4);
end

% change to lower case
monthstr = lower(monthstr);

% handle ambiguous case in French
if strcmp(monthstr,'juil') == 1
    monthstr = 'jul';
elseif strcmp(monthstr,'juin') == 1
    monthstr = 'jun';
else
    monthstr = monthstr(1:3);
end

switch monthstr
    case 'jan'
        monthnum = 1;
    case {'feb','fev'}
        monthnum = 2;
    case 'mar'
        monthnum = 3;
    case {'apr','avr'}
        monthnum = 4;
    case {'may','mai'}
        monthnum = 5;
    case 'jun'
        monthnum = 6;
    case 'jul'
        monthnum = 7;
    case {'aug','auo'}
        monthnum = 8;
    case 'sep'
        monthnum = 9;
    case 'oct'
        monthnum = 10;
    case 'nov'
        monthnum = 11;
    case 'dec'
        monthnum = 12;
    otherwise
        error(sprintf('Unrecognized monthnum %s',monthstr));
        monthnum = 0;
end
return


