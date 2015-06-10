function isnewer=fnewer(fname1,fname2)
%function isnewer=fnewer(fname1,fname2)
% check to see if fname2 is newer than fname1

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

if fexist(fname1) == 1 && fexist(fname2) == 1
    D1=dir(fname1);
    D2=dir(fname2);
    
    if D2.datenum > D1.datenum
        isnewer = 1;
    else
        isnewer = 0;
    end
else
    isnewer = 0;
end
return
