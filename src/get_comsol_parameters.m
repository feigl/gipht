function [pnames32, pnames, values, dims, descrs] = get_comsol_parameters(mphname)
%function [pnames32, pnames, values, dims, descrs] = get_comsol_parameters(mphname)
% get names of the parameters in a COMSOL mph file named mphname

verbose = 1;
if verbose == 1;
    fprintf(1,'Loading %s\n',mphname);
end

%model = svartsengi2DaxiV09a;

model = mphload(mphname);

param_names = model.param.varnames;

expr = mphgetexpressions(model.param);
[nvals,ncols] = size(expr);

%nvals = numel(param_names);
pnames = cell(nvals,1);
pnames32 = cell(nvals,1);
descrs = cell(nvals,1);
values = nan(nvals,1);
dims   = cell(nvals,1);

for ip = 1:nvals
    pname = char(param_names(ip));
    pnames{ip} = sprintf('%s',pname);
    descr = char(model.param.descr(param_names(ip)));
    descrs{ip} = sprintf('%s',char(descr));
    
    % look for numerical value
    %value = char(model.param.get(param_names(ip)));
    %values(ip) = sscanf(value,'%f');
    %value = model.param.get(param_names(ip))
    value = char(expr(ip,2));
    % look for square bracket
    idim = strfind(value,'[');
    if idim > 0
        values(ip) = str2double(value(1:idim-1));
    else
        values(ip) = str2double(value);
    end

        
%     nbad = numel(strfind(value,'*')) ...
%         + numel(strfind(value,'/')) ...
%         + numel(strfind(value,'('));
%     if numel(idim) > 0  && nbad == 0
%         dims{ip} = sprintf('%s',value(idim:end));
%     else
%         dims{ip} = sprintf(' ');
%     end

    % translate special characters
    if idim > 0
    dimstr = value(idim:end);
    dimstr = strrep(dimstr,'*','.');
%     dimstr = strrep(dimstr,'/','DIV');
%     dimstr = strrep(dimstr,'(',' ');
    else
        if isfinite(values(ip)) == 1
            dimstr = ' ';
        else
            dimstr = ' ';
        end
    end
    dims{ip} = dimstr;
    
    
    
    if strcmp(dims{ip},'[1]') == 1
        units_without_brackets = sprintf('dimless');
    else
        units_without_brackets = strrep(strrep(dims{ip},'[',''),']','');
    end
    pname32 = sprintf('%-32s',sprintf('CS %s in %s %s'...
        ,char(pnames{ip}),units_without_brackets...
        ,char(descrs{ip})));
    
    % eliminate special characters
    pname32 = strrep(pname32,'^',' ');
    pname32 = strrep(pname32,'=',' ');
    pname32 = strrep(pname32,'(',' ');
    pname32 = strrep(pname32,')',' ');
    pname32 = strrep(pname32,'in dimless','   dimless');
    pnames32{ip} = strrep(pname32(1:32),' ','_');
end


if verbose == 1
    for ip = 1:nvals
        fprintf(1,'    %-8s %12.4e %s\n',pnames32{ip},values(ip),descrs{ip});
    end
end

return




