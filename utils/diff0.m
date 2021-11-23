function dy = diff0(y)
%function dy = diff0(y)
% return diff with respect to first element
% such that first element of vector is zero and
% dimension of output vector is the same as that of input
% 2021/10/18 Kurt Feigl
dy = zeros(size(y));
dy(2:end) = diff(y);  
return
end

