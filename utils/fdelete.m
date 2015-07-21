function ierr=fdelete(fname)
% % delete a file named fname 

% find out if file exists
if exist(fname,'file')  == 2
    ok = 1;
else
    ok = 0;
end

if ok == 0
    return; % no need to delete if it file does not exist
else
    cmd1=sprintf('/bin/rm -fv %s\n',fname);
    
    %fprintf(1,'Running unix command line:\n%s\n',cmd1);
    [unixstat,unixout] = unix(cmd1);
    
    if unixstat == 0
        %fprintf(1,'Successfully deleted file named %s\n',fname);
        ierr = 0;
    else
        fprintf(1,'ERROR: Could not delete file named %s\n',fname);
        error(unixout);
        ierr = 0;
    end
end
return

