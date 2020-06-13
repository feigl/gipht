function Tparams = get_comsol_parameters2(model,verbose)
%function Tparams = get_comsol_parameters(model,verbose)
% given a comsol model object, return a table of parameters
% 20200507 - supercedes version in gipht
% return a table

narginchk(0,2)
nargoutchk(0,1);

if nargin < 1
    help(sprintf('%s',mfilename));
    return
end
if nargin < 2
    verbose = 1;
end

param_names = model.param.varnames;

expr = mphgetexpressions(model.param);
[mParameters,ncols] = size(expr)

parameterNames = model.param.varnames;

%nvals = numel(param_names);
idNumbers = zeros(mParameters,1);
pnames = cell(mParameters,1);
pnames32 = cell(mParameters,1);
descrs = cell(mParameters,1);
values = nan(mParameters,1);
units   = cell(mParameters,1);


if numel(parameterNames) ~= mParameters
    error('miscount');
end

if verbose == 1
    for i=1:numel(parameterNames)
        try
            value1 = model.param.evaluate(sprintf('%s',parameterNames(i)));
        catch Merror
            Merror
            value1 = nan;
        end
        fprintf(1,'%3d %-32s %10.4G\n',i,parameterNames(i),value1)
%         try
%             unit1 = model.param.evaluateUnit((sprintf('%s',parameterNames(i))));
%         catch Merror
%             Merror
%             unit1 = nan;
%         end
%         
%         fprintf(1,'%3d %-32s %10.4G [%s]\n',i,parameterNames(i),value1,unit1)
     end
end

% get numerical values
for i=1:mParameters
   values(i) = model.param.evaluate(sprintf('%s',parameterNames(i)));
end

% deal with strings in cells
for i = 1:mParameters
    idNumbers(i) = i;
    %  names
    pname1 = char(param_names(i));
    pnames{i} = sprintf('%s',pname1);
    
    % description
    descr1 = char(model.param.descr(param_names(i)));
    descrs{i} = sprintf('%s',char(descr1));
    
%     % units
%     try
%         unit1 = char(model.param.evaluateUnit((sprintf('%s',parameterNames(i)))));
%     catch Mexcept
%         Mexcept
%         unit1 = '';
%     end
    unit1 = 'NaN';

    % handle dimensionless quantities
    if strcmp(unit1,'1') == 1 || numel(unit1) == 0
        unit1 = 'dimensionless';
    end
    unit1 = strrep(unit1,'[','');
    unit1 = strrep(unit1,']','');
    units{i} = unit1;

    % make long name
    pname32 = sprintf('%-32s',sprintf('CS%03d_%s_%s_%s',i,pname1,unit1,descr1));
    
    % eliminate special characters
    pname32 = strrep(pname32,'^',' ');
    pname32 = strrep(pname32,'=',' ');
    pname32 = strrep(pname32,'(',' ');
    pname32 = strrep(pname32,')',' ');
    %pname32 = strrep(pname32,'in dimless','dimless');
    pname32 = strrep(pname32,' ','_');
    pname32 = matlab.lang.makeValidName(pname32);

    % truncate and store
    pnames32{i} = pname32(1:32);
end

Tparams = table(idNumbers,'VariableNames',{'idnumber'});
Tparams = [Tparams,table(pnames32,'VariableNames',{'name32'})];
Tparams = [Tparams,table(pnames,'VariableNames',{'name'})];
Tparams = [Tparams,table(values,'VariableNames',{'value'})];
Tparams = [Tparams,table(units,'VariableNames',{'units'})];
Tparams = [Tparams,table(descrs,'VariableNames',{'description'})];

if verbose == 1
    Tparams
end

return
end




