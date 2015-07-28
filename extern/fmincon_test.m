function x = fmincon_test

%     To optimize for specific values of a1 and a2, first assign the values
%     to these two parameters. Then create two one-argument anonymous
%     functions that capture the values of a1 and a2, and call myfun and
%     mycon with two arguments. Finally, pass these anonymous functions to
%     fmincon:

a1 = 2; 
% a2 = 1.5; % define parameters first
%options = optimset('Algorithm','interior-point','Display','Iter'); % run interior-point algorithm
options = optimset('Algorithm','sqp','Display','Iter'); % run interior-point algorithm
%x = fmincon(@(x) myfun(x,a1),[1;2],[],[],[],[],[],[],@(x) mycon(x,a2),options)
lb = [1.4 0.9]
ub = [2.1 1.5]

fittingfun = str2func('abs')

x1 = lb(1):0.1:ub(1);
x2 = lb(2):0.1:ub(2);
for j = 1:numel(x1)
    for i = 1:numel(x2)        
        %F(i,j) = myfun([x1(j),x2(i)],a1);
        F(i,j) = myfun2([x1(j),x2(i)],a1,fittingfun);
    end
end
size(F)
close all;
figure;
imagesc(x1,x2,F);
axis xy; 
hold on; 
colorbar;

p0 = [1.70, 1.30]; % bad initial estimate   
%p0 = [1.40, 0.90]; % perfect initial estimate
%p1 = fmincon(@(x) myfun(x,a1),p0,[],[],[],[],lb,ub,[],options)
p1 = fmincon(@(x) myfun2(x,a1,fittingfun),p0,[],[],[],[],lb,ub,[],options)

plot(p0(1),p0(2),'ok');
plot(p1(1),p1(2),'*k');
legend('initial','final');

return


function f = myfun(x,a1)
f = x(1)^2 + a1*x(2)^2;
return

function f = myfun2(x,a1,fittingfun)
f = x(1)^2 + fittingfun(a1)*x(2)^2;
return

% function [c,ceq] = mycon(x,a2)
% c = a2/x(1) - x(2);
% ceq = [];
% return;
