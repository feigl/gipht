function [xvalley,yedges] = fitvalley(xvals,yvals)
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

p0(1) = nanmean(xvalley);
p0(2) = (quantile(xvalley,0.975) - quantile(xvalley,0.25))/2.;
lb(1) = nanmin(xvalley); 
ub(1) = nanmax(yvalley);
lb(2) = 0.;
ub(2) = (nanmax(yvalley) - nanmin(xvalley))/2.; 

params = p0;

%[PST.p1,cost1,exitflag,output] = fmincon(@(params) obj_function(params,DST),PST.p0,A,b,Aeq,beq,PST.lb,PST.ub);

% Assign values to the parameters and define a function handle f to an anonymous function
fanon = @(params, xvalley, yvalley) obj_function(params, xvalley, yvalley)

% X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
%     the default optimization parameters replaced by values in OPTIONS, an
%     argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
    
[p1,cost1,exitflag,output] = fmincon(fanon,p0,A,b,Aeq,beq,lb,ub);
cost1
mse1=sqrt(cost1)

xvfit = linspace(nanmin(xvalley),nanmax(yvalley),100);
yvfit = fun_gaussian_bowl(p1(1),p1(2),xvfit);


whos
figure;
hold on;
plot(xvals,yvals,'k+');
plot(xvalley,yvalley,'r-');


title('bins');

else
    error('miscount');
end

return
end

function cost = obj_function(params,x,y)

mu = params(1);
sigma = params(2);

% residual
res = colvec(fun_gaussian_bowl(mu,sigma,x) - y);

cost = res' * res;

return
end

function y = fun_gaussian_bowl(x,mu,sigma)
% return an upside-down Gaussian with mean mu and standard deviation, evaluate at x
% 2020/05/26 Kurt Feigl

y = normpdf(0,mu,sigma)-normpdf(x,mu,sigma);

return
end



