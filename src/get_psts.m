function [P,pnames] = get_psts(fnames)
%function [P,pnames] = get_psts(fnames)
% given a list of PST files, return a structure of parameters
% input: fnames == cell variable containing file names of PST files
%                  e.g. fnames = flist('../gipht21/PSELECT7_20111001/p*/PST.OUT');
% output: 
%        P       == structure of PST structures
%        pnames  == character array with unique list of parameter names
% 2011-OCT-04 Kurt Feigl

k= 0;
npairs=0;
%P = struct([]);
nfiles = numel(fnames);
for i=1:nfiles
    if fexist(fnames{i})
        k=k+1;
        % read one parameter file
        P1 = read_pst(fnames{i});
        % read list of paramter names
        if k==1
            pnames = P1.names;
        end
        % check parameter names
        for j=1:numel(pnames)
            if strcmp(char(pnames{j}),char(P1.names{j})) == 0
                fprintf(1,'WARNING parameter %d does not match. EXPECTED:%s\n FOUND : %s\nEXPECTED:%s\n'...
                    ,char(pnames{j}),char(P1.names{j}));
            end
        end
        P(k)=P1;
    else
        warning(sprintf('Could not open file named %s\n',char(fnames{i})));
    end
end
% figure out the names of parameters
pnames={P.names};
pnames=pnames(1);
pnames=pnames{1};

npairs = numel(P);
if npairs < 1
    error(sprintf('npairs = %d\n',npairs));
end
return

