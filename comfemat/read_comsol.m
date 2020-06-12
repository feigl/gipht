%% read an MPH file output by COMSOL
% 20200420 Kurt Feigl

%% initialize
close all
clear variables
verbose = 0
if verbose == 1
    echo on;
    diary(sprintf('%s_%s.log',mfilename,datestr(now,30)));
else
    echo off;
end
    

%% Comsol server must be started
status = which('mphload');
if contains(status,'mphload') ~= 1
    mphstart
    import com.comsol.model.util.*
else
    fprintf(1,'Matlab Livelink server is available.\n');
    fprintf(1,'Try the following:\n');
    fprintf(1,'  help mphload\n');
    fprintf(1,'  mphnavigator\n');
    fprintf(1,'  mphgetproperties\n'); 
end


%% set path for GIPhT utlities
home = getenv('HOME');
if numel(home) > 0
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'GraphTreeTA')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'extern')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'utils')),'-begin');
    addpath(genpath(strcat(home,filesep,'gipht',filesep,'src')),'-begin');
else
    addpath(genpath('~feigl/gipht/GraphTreeTA'),'-begin');
    addpath(genpath('~feigl/gipht/extern'),'-begin');
    addpath(genpath('~feigl/gipht/utils'),'-begin');
    addpath(genpath('~feigl/gipht/src'),'-begin');
end

%% load the MPH file
fileNameMPH = '01/LdM_3DFSI_spheroid_Tvisco_OPT2019DEC02WORKS_01.mpho' 
%fileNameMPH = '01/LdM_3DFSI_spheroid_Tvisco_OPT2019DEC02WORKS_01withExport.mph' 
model=mphload(fileNameMPH);
info0 = mphsolutioninfo(model)

%mphgetproperties(model)

%% get a table
% we have to know its tag!
Tstruct = mphtable(model,'tbl4');
% get the data
ArrayData = Tstruct.data;
[nrows,ncols] = size(ArrayData)
% get the headers and make them into legitimate variable names 
varnames = Tstruct.headers;
varnames = matlab.lang.makeValidName(varnames);
% make a table
Ttable = array2table(ArrayData);
% assign names to variables
Ttable.Properties.VariableNames =varnames;su

return

%% get information about solutions
nSolutions = numel(info0.solutions)
SolTags = cell(nSolutions,1);
for iSolution = 1:nSolutions
    SolTags{iSolution} = char(info0.solutions{iSolution});
end

for iSolution = 1:nSolutions
    SolInfos{iSolution} = mphsolinfo(model,'soltag',SolTags{iSolution});
end

% look at these carefully
SolInfos{1:nSolutions}

% in our case, the 4th solution contains information about the GPS data
SolInfos{4}.solpar
% 
ttGPS = SolInfos{4}.solvals(1:4:end); % time tags
u1GPS = SolInfos{4}.solvals(2:4:end); % vector component
u2GPS = SolInfos{4}.solvals(3:4:end); % vector component
u3GPS = SolInfos{4}.solvals(4:4:end); % vector component

% look at the parameters for the sweep
SolInfos{4}.paramsweepnames
SolInfos{4}.paramsweepvals
for i=1:numel(SolInfos{4}.paramsweepnames)
    fieldName1 = SolInfos{4}.paramsweepnames{i};
    S4params.(fieldName1) = SolInfos{4}.paramsweepvals(i);
end


%% get Comsol evaluation solution_epochs in seconds
for i=1:nSolutions
    info1 = SolInfos{i};
    tSol1 = info1.solvals;
    t1 = min(tSol1);
    t2 = max(tSol1);
    if (t2 - t1) > 0. 
        fprintf(1,'solution number %in2d tag = "%s" label = "%s"\n' ...
        ,i,info1.soltag,info1.label);
    end
end

return

%% obtain nodal coordinates first:
nodestruct = mphxmeshinfo(model,'soltag','sol1');
coords = nodestruct.nodes.coords;
%max(max(coords));
[ndim,ncols] = size(coords);
x_node = coords(1,:); % x coordinate  in meters, w.r.t. xcen
y_node = coords(2,:); % y coordinate  in meters, w.r.t. ycen
z_node = coords(3,:); % z component in meters, positive upwards, w.r.t. 0
%obtain nodal coordinates first:
nodestruct = mphxmeshinfo(model,'soltag','sol1');
coords = nodestruct.nodes.coords;
%max(max(coords));
[ndim,ncols] = size(coords);
x_node = coords(1,:); % x coordinate  in meters, w.r.t. xcen
y_node = coords(2,:); % y coordinate  in meters, w.r.t. ycen
z_node = coords(3,:); % z component in meters, positive upwards, w.r.t. 0
figure
plot(x_node,y_node,'r+');
figure
nanmax(z_node)
nanmin(z_node)
itop=find(abs(z_node-nanmax(z_node))< eps);
plot(x_node(itop),y_node(itop),'r+');
info1
info2