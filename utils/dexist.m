function ok=dexist(dirname)
% check to see if a directory named dirname exists
% function ok=dexist(dirname)

% Use built-in function
if exist(dirname,'dir')  == 7
    ok = 1;
else
    ok = 0;
end
return
