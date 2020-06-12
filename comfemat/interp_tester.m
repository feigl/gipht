%Example 2: Evaluation including several solution
% from mphdoc

modelName = 'model_tutorial_llmatlab';
fileNameMPH = strcat(modelName,'.mph');
if exist(fileNameMPH,'file') == 0
    % Load model_tutorial_llmatlab.mph:
    model = mphopen();
    % Create a study combining a parametric sweep and a transient
    % study step:
    std = model.study.create('std');
    param = std.feature.create('param','Parametric');
    time = std.feature.create('time','Transient');
    % Set the time stepping and the parametric sweep parameters:
    time.set('tlist', 'range(0,1,25)');
    param.setIndex('pname','power',0).setIndex('plistarr','30 60 90', 0);
    % Run the study:
    std.run;
    
    mphsave(model,fileNameMPH);  
end

model = load_comsol_mph_file(fileNameMPH,1);

% % Evaluate the temperature at every time step computed with power set
% % to 30:
coord = [0 0 1e-2;0 0 1e-2;0 1e-2 1e-2];
T = mphinterp(model,'T','coord',coord,'dataset','dset2')
% 
% % Evaluate the temperature at the fifth time step:
% T = mphinterp(model,'T','coord',coord,'dataset','dset2','solnum',5)
% 
% % Evaluate the temperature at 10.5 sec:
% T = mphinterp(model,'T','coord',coord,'dataset','dset2','t',10.5)
% 
% % Evaluate the temperature at every time step computed with power set
% % to 90:
% T = mphinterp(model,'T','coord',coord,'dataset','dset2',...
%     'outersolnum',3)

% get time steps for solution
info1=mphsolinfo(model,'soltag','sol1')
tsol=info1.solvals
[nTimes,nDummy] = size(tsol)


%obtain nodal coordinates first:
nodestruct = mphxmeshinfo(model,'soltag','sol1');
coords = nodestruct.nodes.coords;
%max(max(coords));
[nDim,nPoints] = size(coords)

Tall = mphinterp(model,{'T'},'coord',coords,'t',tsol); %
[nTrows,nTcols] = size(Tall)

if nTrows ~= nTimes || nTcols ~= nPoints
    error('Dimensions are incorrect');
end

% make a figure;
figure;
plot(tsol,mean(Tall,2),'k+');
xlabel('t');
ylabel('T');


