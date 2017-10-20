function ft = timefun_riser(tepochs, treference1, treference2)
% return value of rising time function f(t)
% 20170915 Kurt Feigl

if treference1 < treference2
   ft = heaviside1(tepochs - treference1) .* (tepochs - treference1) ...
      - heaviside1(tepochs - treference2) .* (tepochs - treference2);
else
    treference1
    treference2
    error('Reference times out of order.\n');
end
return
end
