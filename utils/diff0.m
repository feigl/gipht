function dy = diff0(y)
%function dy = diff0(y)
% return diff with respect to first element
% such that first element of vector is zero and
% dimension of output vector is the same as that of input
% 2021/10/18 Kurt Feigl
% 2024/07/01 handle datetime class
if isdatetime(y)
    dy=repmat(seconds(0),size(y));
    dy(2:end) = diff(y);
else
    dy = zeros(size(y));
    dy(2:end) = diff(y);
    return
end

