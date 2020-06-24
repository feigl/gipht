function [xvalley,yedges] = fitvalley2(xvals,yvals)
% given a set of values, find the valley function (bottom envelope)
% 20200526 Kurt Feigl

narginchk(2,2);
nargoutchk(2,2);

if numel(xvals) == numel(yvals)

nvals = numel(xvals);

% [N,edges,bin] = histcounts(___) also returns an index array, bin, using any of the previous syntaxes. bin is an array of the same size as X whose
% elements are the bin indices for the corresponding elements in X. The number of elements in the kth bin is nnz(bin==k), which is the same as N(k).
[ncounts,xedges,ibin] = histcounts(xvals);

% N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1). The 
%     last bin will also include the right edge such that N(end) will count
%     X(i) if EDGES(end-1) <= X(i) <= EDGES(end).
 
nbins = numel(xedges)-1;
yedges = nan(size(xedges));
for i=1:nbins
    ibin1 = find(ibin == i);
    yedges(i) = nanmin(yvals(ibin1));
end

% find values at centers of bins
xvalley = xedges(1:end-1) + diff(xedges)/2.;
yvalley = nan(size(xvalley));
for i=1:nbins
    yvalley(i) = nanmin([yedges(i),yedges(i+1)]);
end
xvalley=colvec(xvalley);
yvalley=colvec(yvalley);
    
%test obj func2
% cost0=obj_function(PST.p0,DST)
% mse0=sqrt(cost0)

A=[];
b=[];
Aeq=[];
beq=[];
%objective function calls fitting function calc_dis
%fun=@(params) obj_function(params,DST)
%optoptions = optimset('Algorithm','sqp','Display','Iter'); 

PST.p0(1) = nanmean(xvalley);
PST.p0(2) = (quantile(xvalley,0.975) - quantile(xvalley,0.25))/2.;
PST.lb(1) = nanmin(xvalley); 
PST.ub(1) = nanmax(yvalley);
PST.lb(2) = 0.;
PST.ub(2) = (nanmax(xvalley) - nanmin(xvalley))/2.; 

PST.fitfun = @fun_gaussian_bowl;

DST.x   = xvalley;
DST.obs = yvalley;
DST.sigma = ones(size(DST.obs));

TST = [];

params = PST.p0;

opt_options = optimset('MaxFunEvals',1e6); 

[PST.p1,cost1,exitflag,Output] = simulannealbnd(@(params) obj_function2(params,PST,DST,TST),PST.p0...
    ,PST.lb,PST.ub...
    ,opt_options);
cost1
mse1=sqrt(cost1)

% simulated values
SST.obs = colvec(yvalley);
SST.sigma = ones(size(SST.obs));

SST.x   = linspace(nanmin(xvalley),nanmax(xvalley),100);
SST.mod = fun_gaussian_bowl(PST.p1,PST,SST,TST);


whos
figure;
hold on;
plot(xvals,yvals,'k+');
plot(xvalley,yvalley,'b-');
plot(SST.x,SST.mod,'r-');


title('bins');

else
    error('miscount');
end

return
end


function y = fun_gaussian_bowl(params,PST,DST,TST)
% return an upside-down Gaussian with mean mu and standard deviation, evaluate at x
% 2020/05/26 Kurt Feigl

mu = params(1);
sigma = params(2);

y = normpdf(0,mu,sigma)-normpdf(DST.x,mu,sigma);

return
end

function misfit=obj_function2(params, PST, DST, TST)
PST.p0=params;
% get a handle on fitting function
hfun = PST.fitfun;
if isa(hfun,'function_handle') == false
    fprintf(1,'Making function handle\n');
    hfun = @PST.fitfun;
end
% evaluate fitting function with current values of parameters
DST.mod = hfun(params,PST,DST,TST);
DST.resid=DST.obs-DST.mod;
DST.nresid=DST.resid./DST.sigma;

% L1 norm
%misfit=sum(abs(DST.nresid))/numel(DST.nresid);
% L2 norm
misfit=DST.nresid'* DST.nresid;
end

