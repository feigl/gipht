function V  = range2vector(R,S,M,threshm)
%function V  = range2vector(R,S,M,threshm)
% given range, find vector
%    inputs:
%       r = scalar range change in meters
%       s = unit vector pointing from target to sensor
%       m = displacment vector [east, north, up] from model calculation
%       threshm = threshold in meters
%
% 2012-JUN-25 Kurt Feigl
% 2012-OCT-24 new approach:
%           assume V  is parallel to M and scale by ratio of ranges
% 2012-OCT-25 corrected with help from Helene Le Mevel
% 2012-OCT-28 apply threshold to avoid dividing by zero
%
% Example test case:
%
% The following input
%
% R = [1,2,3];
% S = [1,0,0;0,1,0;0,0,1]
% M = 2* S;
% V = range2vector(R,S,M)
%
% should return:
%
% V =
%     -1     0     0
%      0    -2     0
%      0     0    -3
%
nargchk(3,4,nargin);

% default value for threshhold
if exist('threshm','var') == 0
    threshm = 3.0e-3;
end

ndat = numel(R);
r = R;

[nrows,mcols] = size(M);

if numel(S) == 3*ndat && numel(M)== 3*ndat
    s = reshape(S,3,ndat);
    m = reshape(M,3,ndat);
    V  = nan(3,ndat);
    mhat = nan(3,ndat);
    for i=1:ndat
        
        % modeled range in meters
        rm = -1.0 * [m(1,i),m(2,i),m(3,i)] * [s(1,i),s(2,i),s(3,i)]';
        
        % avoid explosion
        if abs(rm) >= threshm
            %range ratio of observed to modeled
            ratio = r(i) / rm;
        else
            ratio = NaN;
        end
        
        
        % east component
        V(1,i) = ratio * m(1,i);
        % north component
        V(2,i) = ratio * m(2,i);
        % upward component
        V(3,i) = ratio * m(3,i);
    end
else
    error('Argument count mismatch');
end
return

end

