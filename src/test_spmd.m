% test spmd function

% parpool(2)
% spmd
%   % build magic squares in parallel
%   q = magic(labindex + 2);
% end
% for ii=1:length(q)
%   % plot each magic square
%   figure, imagesc(q{ii});
% end
% delete(gcp)

parpool(2)
spmd
    fid=fopen(sprintf('test%1d.dat',labindex),'w+');
    fprintf(fid,'labindex is %d\n',labindex);
    fclose(fid);
end
delete(gcp)
