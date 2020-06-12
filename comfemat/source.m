function q=source(x,y)
% vector inputs and outputs

for i=1:length(x)
    q(i) = calcsource(x(i), y(i));
end


function q=calcsource(x,y)
% scalar inputs and outputs

r = sqrt(x^2+y^2);
if r<0.5
    q = 1000;
else
    q = 500;
end
