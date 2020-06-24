function [xvalley,yvalley] = fitvalley3(xvals,yvals)
% given a set of values, find the valley function (bottom envelope)
% 20200526 Kurt Feigl

narginchk(2,2);
nargoutchk(2,2);

if numel(xvals) == numel(yvals)

nvals = numel(xvals);

% find convex hull
dt = delaunayTriangulation(xvals,yvals);
k1 = convexHull(dt); 

% % % only want the bottom half
% kmin = find(abs(xvals-nanmin(xvals)) <= eps)
% kmax = find(abs(xvals-nanmax(xvals)) <= eps)

%k2 = intersect(find(yvals<yvals(kmin)), find(yvals<yvals(kmax))); 
%k2 = find(yvals <= max([yvals(kmin), yvals(kmax)])); 

%k2 = find(yvals <= nanmean(yvals)); 
k2 = find(yvals <= nanmedian(yvals)); 

k = intersect(k1,k2);
%k = unique(k);

if numel(k) < 2
    warning('Too few points');
    return
end

xvalley = xvals(k);
yvalley = yvals(k);

% sort to make clean plot
[xvalley, ksort] = sort(xvalley);
yvalley = yvalley(ksort);

ndata = numel(yvalley)
for i=1:ndata
    fprintf(1,'%10.4f %10.4f\n',xvalley(i),yvalley(i));
end


A=[];
b=[];
Aeq=[];
beq=[];

TST = [];

%objective function calls fitting function calc_dis
%fun=@(params) obj_function(params,DST)
%optoptions = optimset('Algorithm','sqp','Display','Iter'); 

%PST.p0(1) = nanmean(xvalley);
imin = find(abs(yvals-nanmin(yvals)) <= eps)
xvals(imin)
yvals(imin)

PST.p0(1) = xvals(imin);
PST.p0(2) = (quantile(xvalley,0.975) - quantile(xvalley,0.025))/2.;
% PST.lb(1) = nanmin(xvalley); 
% PST.ub(1) = nanmax(xvalley);
PST.lb(1) = xvals(imin); 
PST.ub(1) = xvals(imin);
PST.lb(2) = (nanmax(xvalley) - nanmin(xvalley))/5.; 
PST.ub(2) = (nanmax(xvalley) - nanmin(xvalley))/2.; 

PST.fitfun = @fun_gaussian_bowl;

DST.x   = xvalley;
DST.obs = yvalley;
DST.sigma = ones(size(DST.obs));
DST.mod = fun_gaussian_bowl(PST.p0,PST,DST,TST);


params = PST.p0;

opt_options = optimset('MaxFunEvals',1e6); 

[PST.p1,cost1,exitflag,Output] = simulannealbnd(@(params) obj_function2(params,PST,DST,TST),PST.p0...
    ,PST.lb,PST.ub...
    ,opt_options);
PST.p1
cost1
mse1=sqrt(cost1)
Output

% simulated values
SST.obs = colvec(yvalley);
SST.sigma = ones(size(SST.obs));

SST.x   = linspace(nanmin(xvalley),nanmax(xvalley),100);
SST.mod = fun_gaussian_bowl(PST.p1,PST,SST,TST);

figure;
hold on;
plot(xvals,yvals,'k+');
plot(xvals(k1),yvals(k1),'g-');
plot(xvalley,yvalley,'b-','LineWidth',3);
plot(DST.x,DST.mod,'c-');
plot(SST.x,SST.mod,'r-');
plot(xvals(imin),yvals(imin),'ro','MarkerSize',20,'LineWidth',3); %'MarkerFaceColor','r');


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

