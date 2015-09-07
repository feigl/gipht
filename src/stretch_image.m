function image_out = stretch_image(image_in,nlevels)
% stretch an image based on quantiles
% return indices
%  determine cut points from observed values 
nlevels
cuts = quantile(colvec(image_in),linspace(0.0,1.0,nlevels));
figure;
plot(cuts,'r-');
xlabel('index');
ylabel('cut');


[nrows,ncols] = size(image_in);

image_out = ones(nrows,ncols,'uint8');
%
% % slow, but works
for j=1:ncols
    for i=1:nrows
        qtemp = image_in(i,j);
        if isfinite(qtemp) == 1
            for k=1:nlevels-1
                if qtemp >= cuts(k) && qtemp < cuts(k+1)
                    image_out(i,j) = k;
                    %break;
                end
            end
        end
    end
end
return

%end

