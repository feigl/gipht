function ft = time_function_heaviside(tepochs, tquake)
% return value of time function f(t)
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tquake  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch
%                    Heaviside step function
%          ft(t)    = 1 if t >= tquake
%          ft(t)    = 0 if t <  tquake
%
% 2011-OCT-17 Kurt Feigl
[me, ncols] = size(tepochs);
if ncols == 1
    ft = zeros(me,1);
    
    if tquake >= nanmin(tepochs) && tquake <= nanmax(tepochs)
        %fprintf(1,'Using step function that turns on at epoch %f\n',tquake);
        itime=find(tepochs >= tquake);
        ft(itime) = 1.0;
    end
else   
    ncols
    error('Dimension problem');
end
return
end

