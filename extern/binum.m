function [val,spl,sgf] = binum(str,uni)
% Convert a binary prefixed string into numeric values. (ISO/IEC 80000-13)
%
% (c) 2014 Stephen Cobeldick
%
% ### Function ###
%
% Convert a string (containing numeric coefficients with binary prefixes) into
% the equivalent numeric values, and also return the split string parts.
% Identifies both symbol prefixes and full names, e.g. either 'Ki' or 'kibi'.
%
% Syntax:
%  val = binum(str)           % Allow any units in <str>.
%  val = binum(str,uni)       % Define the units in <str>.
%  [val,spl,sgf] = binum(...) % Return the values, the split parts of <str>, & sig-figs.
%
% See also BIPRE SINUM SIPRE STR2NUM NUM2STR MAT2STR SSCANF SPRINTF ROUND60063 ROUND2SF ROUND2DP NUM2WORDS
%
% ### Examples ###
%
% binum('10 Ki')  OR  binum('10.0 kibi')  OR  binum('10240')  OR  binum('1.024e4')
%   ans = 10240
%
% [val,spl] = binum('Memory: 200 mebibyte')
%   val = 209715200
%   spl = {'Memory: ','byte'}
%
% [val,spl,sgf] = binum('From -3.6 MiB to +1.24KiB data transfer.')
%   val = [-3774873.6,1269.8]
%   spl = {'From ','B to ','B data transfer.'}
%   sgf = [2,3]
%
% [val,spl] = binum('100 Pixel','Pixel') % Try it without the second option.
%   val = 100
%   spl = {'','Pixel'}
%
% binum(bipre(9*1024^4))
%   ans = 9895604649984 = 9*1024^4
%
% ### String Format ###
%
% - Any number of coefficients may occur in the string.
% - The coefficients may be any combination of digits, positive or negative,
%   integer or decimal, exponents may be included using E-notation (e/E). 
% - An Inf or NaN value in the string will also be converted to a numeric.
% - The space-character between the coefficient and the prefix is optional.
% - The prefix is optional, either as the binary prefix symbol or name.
% - By default checks first for prefix names, then symbols.
%
% Optional input <uni> controls the prefix/units recognition: if the units may
% contain the prefix characters, then this argument should be specified.
%
% ### Binary Prefix Strings (ISO/IEC 80000-13) ###
%
% Order  |1024^1 |1024^2 |1024^3 |1024^4 |1024^5 |1024^6 |1024^7 |1024^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kibi  | mebi  | gibi  | tebi  | pebi  | exbi  | zebi  | yobi  |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol |  Ki   |  Mi   |  Gi   |  Ti   |  Pi   |  Ei   |  Zi   |  Yi   |
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
% [val,spl,sgf] = binum(str,*uni)

% ### Input Wrangling ###
%
nam = 'kibi|mebi|gibi|tebi|pebi|exbi|zebi|yobi';
sym = 'Ki|Mi|Gi|Ti|Pi|Ei|Zi|Yi';
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
        prc = {'kibi','mebi','gibi','tebi','pebi','exbi','zebi','yobi';...
               'Ki',  'Mi',  'Gi',  'Ti',  'Pi',  'Ei',  'Zi',  'Yi'};
        [~,col] = cellfun(@(s)find(strcmp(s,prc)),tkn(3,idx));
        % Adjust values by coefficients:
        val(idx) = pow2(val(idx),10*col);
    end
    if nargout>2
        sgf = cellfun(@(s)sum(isstrprop(s,'digit')),tkn(1,:));
    end
end
%
end
%----------------------------------------------------------------------END:binum