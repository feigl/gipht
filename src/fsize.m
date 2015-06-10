function isize=fsize(fname)
%function isize=fsize(fname)
% return the size in bytes of file named fname

% DIR List directory.
%     DIR directory_name lists the files in a directory. Pathnames and
%     wildcards may be used.  For example, DIR *.m lists all the M-files
%     in the current directory.
%  
%     D = DIR('directory_name') returns the results in an M-by-1
%     structure with the fields: 
%         name    -- Filename
%         date    -- Modification date
%         bytes   -- Number of bytes allocated to the file
%         isdir   -- 1 if name is a directory and 0 if not
%         datenum -- Modification date as a MATLAB serial date number.
%                    This value is locale-dependent.

if fexist(fname) == 1
    D=dir(fname);
    isize=D.bytes;
else
    isize=0;
end
return
