function h = hist2(d,nbin)
% function h = hist2(d,nbin)
% make a histogram that imitates hist command that used to work 
% in releases of Matlab prior to R2009B
% inputs:
%      d = data (vector)
%      nbin = number of bins (scalar, optional)
% returns:
%      h = graphics handle
%
% Kurt Feigl 2010-MAR-22
% 
% Use this function to work around the (annoying) error messages below
%
% ??? Error using ==> feval
% Undefined function or variable 'datamanager.linkbehavior'.
% 
% Error in ==> hgbehaviorfactory>localCreate at 45
%          ret_h(end+1) = feval(info.constructor);
% 
% Error in ==> hgbehaviorfactory at 31
%     b = localCreate(behavior_name);
% 
% Error in ==> hggetbehavior>localGet at 90
%         b = hgbehaviorfactory(bn);
% 
% Error in ==> hggetbehavior at 72
%     ret_h = localGet(h,behavior_name);
% 
% Error in ==> hist at 116
%       linkBehavior = hggetbehavior(histPatch(k),'Linked');
%  

if nargin <= 2
    nbin = 10;
end
x = linspace(min(d),max(d),nbin);
x = x(1:numel(x)-1) + diff(x)/2;
binned = histc(d,x);
h = bar(x,binned,1.0);
set(h,'FaceColor','blue');
axis tight;
return

