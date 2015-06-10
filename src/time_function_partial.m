function ft = time_function_partial_partial(type, tepochs, tref, y0)
% return value partial derivative of time function f(t)
%
% 2011-OCT-17 Kurt Feigl
%
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tref  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch
%          type == 'heaviside' Heaviside step function
%             ft(t)    = 0 if t <  tref
%             ft(t)    = 1 if t >= tref
%
% Example
% 
% t = 2000:0.1:2010;figure;hold on;
% plot(t,time_function_partial('step',t,2005),'ro-');
% plot(t,time_function_partial('rate',t,2005),'g+-');
% plot(t,time_function_partial('pwl',t,[2005 2007]),'k*-');
% plot(t,time_function_partial('pwl',t,[min(t) 2005 max(t)]),'ms-');
% legend('step','rate','pwl','pwl2');xlabel('t');ylabel('f(t)');


if nargin < 3
    error(sprintf('wrong number of arguments %d. Need at least 3\n',nargin));
end

ft = zeros(size(tepochs));

switch(lower(type))
    case {'rate','secular'}
        %itime=find(isfinite(tepochs));
        ft = tepochs-tref;
    case {'step','heaviside'}
        %itime=find(tepochs >= tref);
        %ft(tepochs >= tref) = 1.0;
        ft = heaviside1(tref);
    case {'pwl'}
        if ~exist('y0','var');
            y0 = 0.0;
        end
        ft = y0 * ones(size(tepochs));
        for i=1:numel(tepochs)
           % accumulate sum
           %ft(i) = y0;
           for j=1:numel(tref)-1
                %ft(i) = ft(i) + heaviside1(tepochs(i)-tref(j));
                 %if tepochs(i) >= tref(j) && tepochs(i) < tref(j+1)
                  if tepochs(i) >= tref(j+1)
                      ft(i) = ft(i) + (tref(j+1)-tref(j));
                  end
                  if tepochs(i) >= tref(j)
                     %ft(i) = ft(i) + (tepochs(i)-tref(j));
                     ft(i) = ft(i) + (tepochs(i)-tref(j));
                 end
            end
        end
    otherwise
        warning(sprintf('undefined type %s',type));
end
return
end

