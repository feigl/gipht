function r=rbarrad(pha)
%function r=rbarrad(pha1)
% return mean resultant length of phase pha in radians
%
% prune NaN
iok = find(isfinite(pha)==1);
pha=pha(iok);

n = numel(pha);
if n > 0
    cbar=sum(cos(pha))/n;
    sbar=sum(sin(pha))/n;
    r=sqrt(cbar^2 + sbar^2);
else
    warning('Number of data points is less than or equal to zero.');
    r = NaN;
end

if r < 0
    warning('Value of r is less than zero.');
    r = 0;
end
if r > 1
    warning('Value of r is greater than unity.');
    r = 1;
end

return
