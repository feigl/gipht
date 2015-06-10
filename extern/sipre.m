function str = sipre(val,sgf,pfx,trz)
% Convert a scalar numeric into an SI prefixed string. (metric/engineering)
%
% (c) 2014 Stephen Cobeldick
%
% ### Function ###
%
% Convert a scalar numeric value into a string. The value is shown in the string
% as a coefficient and an SI unit prefix, optimally chosen for readability. If the
% rounded |val|<10^-24 or |val|>=10^27 then E-notation is used, without a prefix.
%
% Syntax:
%  str = sipre(val)             % Four significant figures and prefix symbol.
%  str = sipre(val,sgf)         % Select significant figures, prefix symbol.
%  str = sipre(val,sgf,pfx)     % Select sig-figs, choose prefix symbol or name.
%  str = sipre(val,sgf,pfx,trz) % Select if decimal trailing zeros are required.
%
% See also SINUM BIPRE BINUM NUM2STR STR2NUM MAT2STR SSCANF SPRINTF ROUND60063 ROUND2SF ROUND2DP NUM2WORDS
%
% ### Examples ###
%
% sipre(10000)  OR  sipre(1e4)
%   ans = '10 k'
% sipre(10000,4,true)
%   ans = '10 kilo'
% sipre(10000,4,false,true)
%   ans = '10.00 k'
%
% ['Power: ',sipre(200*1000^2,2,true),'watt']
%   ans = 'Power: 200 megawatt'
%
% sipre(-5.555e9,2) % Rounds significant figures correctly.
%   ans = '-5.6 G'
%
% sprintf('Clock frequency is %shertz.',sipre(1234567890,5,true))
%   ans = 'Clock frequency is 1.2346 gigahertz.'
%
% sipre(sinum('9 T'))
%   ans = '9 T'
%
% ### SI Prefix Strings ###
%
% Order  |1000^1 |1000^2 |1000^3 |1000^4 |1000^5 |1000^6 |1000^7 |1000^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kilo  | mega  | giga  | tera  | peta  | exa   | zetta | yotta |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   k   |   M   |   G   |   T   |   P   |   E   |   Z   |   Y   |
%
% Order  |1000^-1|1000^-2|1000^-3|1000^-4|1000^-5|1000^-6|1000^-7|1000^-8|
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | milli | micro | nano  | pico  | femto | atto  | zepto | yocto |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   m   |   u   |   n   |   p   |   f   |   a   |   z   |   y   |
%
% ### Input & Output Arguments ###
%
% Inputs (*=default):
%  val = NumericScalar, the value to be converted to string <str>.
%  sgf = NumericScalar, the significant figures in the coefficient, *4.
%  pfx = LogicalScalar, true/false* -> select SI prefix as name/symbol.
%  trz = LogicalScalar, true/false* -> select if decimal trailing zeros are required.
%
% Output:
%  str = Input <val> as a string: coefficient + space character + SI prefix.
%
% str = sipre(val,*sgf,*pfx,*trz)

% ### Input Wrangling ###
%
if nargin<4
    trz = false;
else
    assert(islogical(trz)&&isscalar(trz),'Fourth input <trz> must be a logical scalar.')
end
if nargin<3
    pfx = false;
else
    assert(islogical(pfx)&&isscalar(pfx),'Third input <pfx> must be a logical scalar.')
end
if nargin<2
    sgf = 4;
else
    assert(isnumeric(sgf)&&isscalar(sgf),'Second input <sgf> must be a numeric scalar.')
    sgf = double(uint8(sgf));
end
assert(isnumeric(val)&&isscalar(val)&&isreal(val),'First input <val> must be a real numeric scalar.')
val = double(val);
%
if trz && sgf>1
    fmt = '%#.*g %s';
else
    fmt = '%.*g %s';
end
%
% ### Generate String ###
%
if isfinite(val)
    % Calculate coefficient value:
    xpt = rem(min(9,max(-9,[0;1]+floor(log10(abs(val))/3))),9);
    cof = val.*1000.^-xpt;
    % Round coefficient value:
    ord = 1+floor(log10(abs(cof)));
    if val~=0
        cof = 10.^(ord-sgf).*round(cof.*10.^(sgf-ord));
    end
    % Select prefix symbol/name:
    pfc = {'yocto','zepto','atto','femto','pico','nano','micro','milli',...
        '','kilo', 'mega', 'giga','tera', 'peta','exa', 'zetta','yotta';...
           'y',    'z',    'a',   'f',    'p',   'n',   'u',    'm',...
        '','k',    'M',    'G',   'T',    'P',   'E',   'Z',    'Y'};
    idx = 1+any(abs(cof)==[1000;1]);
    pfs = pfc{2-pfx,9+xpt(idx)};
    % Convert to string (without prefix || digits>whole part):
    if abs(ord(idx))>4 || floor(log10(abs(cof(idx)))-sgf)<-1
        str = sprintf(fmt,sgf,cof(idx),pfs);
    else % (digits<=whole part)
        str = sprintf('%.0f %s',cof(idx),pfs);
    end
else
    str = sprintf('%f ',val);
end
%
end
%----------------------------------------------------------------------END:sipre