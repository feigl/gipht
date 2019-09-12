function ok=fexist(fname)
%function ok=fexist(fname)
% % check to see if a file named fname exists
% returns 0 if fname does not exist
% returns 1 if fname does exist

% home-made version
if numel(fname) > 0
    [fid,message]=fopen(fname,'r');
    if fid > 0
        fclose(fid);
        ok = 1;
    else
        message
        ok = 0;
    end
else
    error(sprintf('File name "fname" not defined.\n'));
    ok = 0;
end

% Use built-in function instead
% % This looks anywhere in the Matlab command search path
% if exist(fname,'file')  == 2
%     ok = 1;
% else
%     ok = 0;
% end
return
