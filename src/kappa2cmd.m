function cmd = kappa2cmd(kappa,n)
%function cmd = kappa2cmd(kappa,n);
% given inputs:
%    kappa: concentration parameter for a von Mises Distribution 
%    n:     sample size
% return outputs:
%    cmd:   circular mean deviation CMD 
%           in radians
%           assuming mean direction of zero
%           defined by Mardia 1972 who denotes it  d0
%
% Kurt Feigl 2008-FEB-10
%
if kappa > 0
   mu = 0.0; % set to zero
   for i=1:1000
      y = randraw('vonmises', [mu, kappa], [1 n]);
      cmd0(i) = circular_mean_deviation(y,0);  % do not remove mean direction
      %cmd0(i) = circular_mean_deviation(y,1);  % remove mean direction
   end
%    figure
%    hist(cmd0)
   %cmd=mean(cmd0);
   cmd=nanmean(cmd0);
else
   warning('undefined kappa')
   cmd = NaN;
end
return

