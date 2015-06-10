function [val,spl,sgf] = sinum(str,uni)
% Convert an SI prefixed string into numeric values. (metric/engineering)
%
% (c) 2014 Stephen Cobeldick
%
% ### Function ###
%
% Convert a string (containing numeric coefficients with SI prefixes) into
% the equivalent numeric values, and also return the split string parts.
% Identifies both symbol prefixes and full names, e.g. either 'k' or 'kilo'.
%
% Syntax:
%  val = sinum(str)           % Allow any units in <str>.
%  val = sinum(str,uni)       % Define the units in <str>.
%  [val,spl,sgf] = binum(...) % Return the values, the split parts of <str>, & sig-figs.
%
% See also SIPRE BINUM BIPRE STR2NUM NUM2STR MAT2STR SSCANF SPRINTF ROUND60063 ROUND2SF ROUND2DP NUM2WORDS
%
% ### Examples ###
%
% sinum('10 k')  OR  sinum('10.0 kilo')  OR  sinum('10000')  OR  sinum('1e4')
%   ans = 10000
%
% [val,spl] = sinum('Power: 200 megawatt')
%   val = 200000000
%   spl = {'Power: ','watt'}
%
% [val,spl,sgf] = sinum('from -3.6 MV to +1.24kV potential difference.')
%   val = [-3600000,1240]
%   spl = {'from ','V to ','V potential difference.'}
%   sgf = [2,3]
%
% [val,spl] = sinum('100 meter','meter') % Try it without the second option.
%   val = 100
%   spl = {'','meter'}
%
% sinum(sipre(9*1000^4))
%   ans = 9000000000000 = 9*1000^4
%
% ### String Format ###
%
% - Any number of coefficients may occur in the string.
% - The coefficients may be any combination of digits, positive or negative,
%   integer or decimal, exponents may be included using E-notation (e/E). 
% - An Inf or NaN value in the string will also be converted to a numeric.
% - The space-character between the coefficient and the prefix is optional.
% - The prefix is optional, either as the SI prefix symbol or name.
% - By default checks first for prefix names, then symbols.
%
% Optional input <uni> controls the prefix/units recognition: if the units may
% contain the prefix characters, then this argument should be specified.
%
% ### SI Prefix Strings ###
%
% Order  |1000^1 |1000^2 |1000^3 |1000^4 |1000^5 |1000^6 |1000^7 |1000^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kilo  | mega  | giga  | tera  | peta  | exa   | zetta | yotta |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol |   k   |   M   |   G   |   T   |   P   |   E   |   Z   |   Y   |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
%
% Order  |1000^-1|1000^-2|1000^-3|1000^-4|1000^-5|1000^-6|1000^-7|1000^-8|
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | milli | micro | nano  | pico  | femto | atto  | zepto | yocto |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol |   m   |   u   |   n   |   p   |   f   |   a   |   z   |   y   |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
%
% ### Input and Output Arguments ###
%
% Inputs (*=default):
%  str = String, with coefficients and prefixes to convert to numeric values.
%  uni = String, to specify the units that are given after the prefix.
%      = LogicalScalar, true/false -> check only for prefix symbol/name.
%      = *[], automagically check for prefix name or symbol, with any units.
%
% Outputs:
%  val = NumericVector, with values calculated from the coefficients and prefixes
%        given in <str>. The size is 1xN, N = number of detected coefficients.
%  spl = CellOfStrings, parts of <str> split by the detected coefficients(+prefixes).
%  sgf = NumericVector, same size as <val>, significant figures of each coefficient.
%
% [val,spl,sgf] = sinum(str,*uni)

% ### Input Wrangling ###
%
nam = 'yocto|zepto|atto|femto|pico|nano|micro|milli|kilo|mega|giga|tera|peta|exa|zetta|yotta';
sym = 'y|z|a|f|p|n|u|m|k|M|G|T|P|E|Z|Y';
sep = '|';
fol = '';
%
% Determine the prefix+unit combination:
if nargin<2||(isnumeric(uni)&&isempty(uni))
    % Name/symbol prefix, any units.
elseif ischar(uni)&&isrow(uni)
    % Units are the given string:
    fol = ['(?=',regexptranslate('escape',uni),')'];
elseif islogical(uni)&&isscalar(uni)
    sep = '';
    if uni % Prefix symbols only.
        nam = '';
    else   % Prefix names only.
        sym = '';
    end
else
    error('Second input <uni> must be a logical scalar, a string, or empty numeric.')
end
assert(ischar(str)&&isrow(str),'First input <str> must be a string.')
%
% ### String Parsing ###
%
% Try to locate a coefficient, possibly with a prefix:
[tkn,spl] = regexp(str,['((+|-)?(NaN|Inf|\d+\.?\d*))(?(1)((e|E)(+|-)?\d+)?) ?',...
    '(',nam,sep,sym,')?',fol],'tokens','split');
%
if isempty(tkn)
    % No coefficient found:
    val = [];
    sgf = [];
else
    % Calculate values from the coefficients:
    tkn = reshape([tkn{:}],3,[]);
    val = cellfun(@(s,x)sscanf([s,x],'%f'),tkn(1,:),tkn(2,:));
    idx = ~cellfun('isempty',tkn(3,:));
    if any(idx)
        % Identify the found prefixes:
        prc = {'yocto','zepto','atto','femto','pico','nano','micro','milli',...
            '','kilo', 'mega', 'giga','tera', 'peta','exa', 'zetta','yotta';...
               'y',    'z',    'a',   'f',    'p',   'n',   'u',    'm',...
            '','k',    'M',    'G',   'T',    'P',   'E',   'Z',    'Y'};
        [~,col] = cellfun(@(s)find(strcmp(s,prc)),tkn(3,idx));
        % Adjust values by coefficients:
        val(idx) = val(idx).*1000.^(col-9);
    end
    if nargout>2
        sgf = cellfun(@(s)sum(isstrprop(s,'digit')),tkn(1,:));
    end
end
%
end
%----------------------------------------------------------------------END:sinum