function ft = timefun_linear(tepochs, treference)
% return value of time function f(t)
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tquake  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch

[me, ncols] = size(tepochs);
if ncols == 1
    ft = zeros(me,1);
    
    if treference >= nanmin(tepochs) && treference <= nanmax(tepochs)
        itime=find(tepochs >= treference);
        ft(itime) = 1.0;
%         fprintf(1,'Using step function that turns on at epoch %f for %d epochs\n',tquake,numel(itime));
    else
        ft = tepochs - treference;
        %fprintf(1,'Linear in time (secular deformation) ft = %10.4f\n',ft);
    end
else
    ncols
    error('Dimension problem');
end
return

