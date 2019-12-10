function y = nicefig(x)
% given a value x, return value rounded to nearest multiple of a power of ten
y=roundn(x,floor(log10(abs(x))));
return


