function str = bipre(val,sgf,pfx,trz)
% Convert a scalar numeric into a binary prefixed string. (ISO/IEC 80000-13)
%
% (c) 2014 Stephen Cobeldick
%
% ### Function ###
%
% Convert a scalar numeric value into a string. The value is shown in the string
% as a coefficient and a binary unit prefix, optimally chosen for readability. If the
% rounded |val|<10^-4 or |val|>=1024^9 then E-notation is used, without a prefix.
%
% Syntax:
%  str = bipre(val)             % Four significant figures and prefix symbol.
%  str = bipre(val,sgf)         % Select significant figures, prefix symbol.
%  str = bipre(val,sgf,pfx)     % Select sig-figs, choose prefix symbol or name.
%  str = bipre(val,sgf,pfx,trz) % Select if decimal trailing zeros are required.
%
% See also BINUM SIPRE SINUM NUM2STR STR2NUM MAT2STR SSCANF SPRINTF ROUND60063 ROUND2SF ROUND2DP NUM2WORDS
%
% ### Examples ###
%
% bipre(10240)  OR  bipre(1.024e4)  OR  bipre(pow2(10,10))  OR  bipre(10*2^10)
%   ans = '10 Ki'
% bipre(10240,4,true)
%   ans = '10 kibi'
% bipre(10240,4,false,true)
%   ans = '10.00 Ki'
%
% ['Memory: ',bipre(200*1024^2,2,true),'byte']
%   ans = 'Memory: 200 mebibyte'
%
% bipre(-5.555e9,2) % Rounds significant figures correctly.
%   ans = '-5.2 Gi'
%
% sprintf('Data saved in %sbytes.',bipre(1234567890,5,true))
%   ans = 'Data saved in 1.1498 gibibytes.'
%
% bipre(binum('9 Ti'))
%   ans = '9 Ti'
%
% ### Binary Prefix Strings (ISO/IEC 80000-13) ###
%
% Order  |1024^1 |1024^2 |1024^3 |1024^4 |1024^5 |1024^6 |1024^7 |1024^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kibi  | mebi  | gibi  | tebi  | pebi  | exbi  | zebi  | yobi  |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|  Ki   |  Mi   |  Gi   |  Ti   |  Pi   |  Ei   |  Zi   |  Yi   |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
%
% ### Input & Output Arguments ###
%
% Inputs (*=default):
%  val = NumericScalar, the value to be converted to string <str>.
%  sgf = NumericScalar, the significant figures in the coefficient, *4.
%  pfx = LogicalScalar, true/false* -> select binary prefix as name/symbol.
%  trz = LogicalScalar, true/false* -> select if decimal trailing zeros are required.
%
% Output:
%  str = Input <val> as a string: coefficient + space character + binary prefix.
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
    xpt = rem(min(9,max(0,[0;1]+floor(log2(abs(val))/10))),9);
    cof = pow2(val,-10*xpt);
    % Round coefficient value:
    ord = 1+floor(log10(abs(cof)));
    if val~=0
        cof = 10.^(ord-sgf).*round(cof.*10.^(sgf-ord));
    end
    % Select prefix symbol/name:
    pfc = {'','kibi','mebi','gibi','tebi','pebi','exbi','zebi','yobi';...
           '','Ki',  'Mi',  'Gi',  'Ti',  'Pi',  'Ei',  'Zi',  'Yi'};
    idx = 1+any(abs(cof)==[1024;1]);
    pfs = pfc{2-pfx,1+xpt(idx)};
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
%----------------------------------------------------------------------END:bipre