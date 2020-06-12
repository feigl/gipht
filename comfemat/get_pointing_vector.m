function scomp = get_pointing_vector(Satellite,Track,icomponent)
%function scomp = get_pointing_vector(Satellite,Track,icomponent)
% given satellite track id, return ith component of unit vector pointing from target on ground to sensor in orbit along line of sight
% suitable for use with COMSOL
% example:
%      Se = get_pointing_vector('ALOS',112,1)
%      Sn = get_pointing_vector('ALOS',112,2)
%      Su = get_pointing_vector('ALOS',112,3)
% Se =
%    -0.5878
% Sn =
%    -0.1686
% Su =
%     0.7912
% 20191110 Kurt Feigl

filename = 'pointing_vectors.xlsx';
T=readtable(filename);
S=table2struct(T,'ToScalar',true);

iok = find(contains(S.Satellite,Satellite));
iok = intersect(iok,find(S.Track == Track));
if numel(iok) == 1
    switch icomponent
        case 1
            scomp = S.Se;
        case 2
            scomp = S.Sn;
        case 3
            scomp = S.Su;
        otherwise
            error(sprintf('Unknown icomponent %d\n',icomponent));
    end
else
    error('problem');
end
return
end

