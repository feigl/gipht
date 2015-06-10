function fnames = flist(match)
% return cell containing file names that match wild card
cmd = sprintf('\\ls -1 %s >! fnames.tmp',match);
%istatus=unix('\ls -1 ../gipht21/*/PST.OUT >! fnames.tmp');
istatus=unix(cmd);
if istatus ~= 0
    istatus
    error('Unix command failed');
end
kount = 1;
k= 0;
fnames = {''};
fid=fopen('fnames.tmp','r');
while kount > 0
    [A,kount] = fscanf(fid,'%s\n',1);
    if kount > 0
        k = k+1;
        fnames{k}=char(A);
    end
end
fclose(fid);
fdelete('fnames.tmp');
return


